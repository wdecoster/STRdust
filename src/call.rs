use bio::io::bed;
use log::{debug, error};
use minimap2::*;
use rayon::prelude::*;
use regex::Regex;
use rust_htslib::{bam, bam::Read};
use std::collections::HashMap;
use std::sync::Mutex;
use url::Url;

use crate::Cli;

pub fn genotype_repeats(args: Cli) {
    debug!("Genotyping STRs in {}", args.bam);
    match (&args.region, &args.region_file) {
        (Some(_region), Some(_region_file)) => {
            error!("ERROR: Specify either a region (-r) or region_file (-R), not both!\n\n");
            panic!();
        }
        (None, None) => {
            error!("ERROR: Specify one of region (-r) or region_file (-R)!\n\n");
            panic!();
        }
        (Some(region), None) => {
            let repeat = crate::repeats::process_region(region)
                .expect("Error when parsing the region string");
            match genotype_repeat(repeat, &args) {
                Ok(output) => {
                    crate::vcf::write_vcf_header(&args.fasta, &args.bam, &args.sample);
                    println!("{output}")
                }
                Err(chrom) => error!("Contig {chrom} not found in bam file"),
            };
        }
        (None, Some(region_file)) => {
            // TODO: check if bed file is okay?
            let mut reader =
                bed::Reader::from_file(region_file).expect("Problem reading bed file!");
            crate::vcf::write_vcf_header(&args.fasta, &args.bam, &args.sample);
            if args.threads > 1 {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(args.threads)
                    .build()
                    .unwrap();
                // chrom_reported and genotypes are vectors that are used by multiple threads to add findings, therefore as a Mutex
                // chrom_reported contains those chromosomes for which an error (absence in the bam) was already reported
                // to avoid reporting the same error multiple times
                let chrom_reported = Mutex::new(Vec::new());
                // genotypes contains the output of the genotyping, a struct instance
                let genotypes = Mutex::new(Vec::new());
                // par_bridge does not guarantee that results are returned in order
                reader.records().par_bridge().for_each(|record| {
                    let rec = record.expect("Error reading bed record.");
                    match genotype_repeat(
                        crate::repeats::RepeatInterval {
                            chrom: rec.chrom().to_string(),
                            start: rec.start().try_into().unwrap(),
                            end: rec.end().try_into().unwrap(),
                        },
                        &args,
                    ) {
                        Ok(output) => {
                            let mut geno = genotypes.lock().unwrap();
                            geno.push(output);
                        }
                        Err(chrom) => {
                            // For now the Err is only used for when a chromosome from the bed file does not appear in the bam file
                            // this error is reported once per chromosome
                            let mut chroms_reported = chrom_reported.lock().unwrap();
                            if !chroms_reported.contains(&chrom) {
                                error!("Contig {chrom} not found in bam file");
                                chroms_reported.push(chrom);
                            }
                        }
                    };
                });
                let mut genotypes_vec = genotypes.lock().unwrap();
                // The final output is sorted by chrom, start and end
                genotypes_vec.sort_unstable();
                for g in &mut *genotypes_vec {
                    println!("{g}");
                }
            } else {
                // When running single threaded things become easier and the tool will require less memory
                // Output is returned in the same order as the bed, and therefore not sorted before writing to stdout
                // chrom_reported contains those chromosomes for which an error (absence in the bam) was already reported
                // to avoid reporting the same error multiple times
                let mut chrom_reported = Vec::new();
                // genotypes contains the output of the genotyping, a struct instance
                for record in reader.records() {
                    let rec = record.expect("Error reading bed record.");
                    match genotype_repeat(
                        crate::repeats::RepeatInterval {
                            chrom: rec.chrom().to_string(),
                            start: rec.start().try_into().unwrap(),
                            end: rec.end().try_into().unwrap(),
                        },
                        &args,
                    ) {
                        Ok(output) => {
                            println!("{output}");
                        }
                        Err(chrom) => {
                            // For now the Err is only used for when a chromosome from the bed file does not appear in the bam file
                            // this error is reported once per chromosome
                            if !chrom_reported.contains(&chrom) {
                                error!("Contig {chrom} not found in bam file");
                                chrom_reported.push(chrom);
                            }
                        }
                    };
                }
            }
        }
    }
}

/// This function genotypes a particular repeat defined by chrom, start and end in the specified bam file
/// All indel cigar operations longer than minlen are considered
/// The bam file is expected to be phased using the HP tag
fn genotype_repeat(
    repeat: crate::repeats::RepeatInterval,
    args: &Cli,
) -> Result<crate::vcf::VCFRecord, String> {
    let flanking = 2000;
    let mut flags = vec![];
    let repeat_ref_seq = match repeat.reference_repeat_sequence(&args.fasta) {
        Some(seq) => seq,
        // Return a missing genotype if the repeat is not found in the fasta file
        None => {
            return Ok(crate::vcf::VCFRecord::missing_genotype(
                &repeat,
                "N",
                ".".to_string(),
            ))
        }
    };

    let repeat_compressed_reference = repeat.make_repeat_compressed_sequence(&args.fasta, flanking);

    let seqs = match get_overlapping_reads(&args.bam, &repeat, args.unphased) {
        Some(seqs) => seqs,
        None => {
            // Return a missing genotype if no (phased) reads overlap the repeat
            return Ok(crate::vcf::VCFRecord::missing_genotype(
                &repeat,
                &repeat_ref_seq,
                0.to_string(),
            ));
        }
    };

    // Create an index for minimap2 alignment to the artificial reference
    let aligner = minimap2::Aligner::builder()
        .map_ont()
        .with_cigar()
        .with_seq(&repeat_compressed_reference)
        .unwrap_or_else(|err| panic!("Unable to build index:\n{err}"));
    let mut consenses: Vec<crate::consensus::Consensus> = vec![];
    let mut all_insertions = vec![]; // only used with `--somatic`
    if args.unphased {
        let mut insertions = vec![];
        // get the sequences of this phase (or all if unphased)
        let seq = seqs.get(&0).unwrap();
        debug!("{repeat}: Unphased: Aliging {} reads", seq.len());
        // align the reads to the new repeat-compressed reference
        for s in seq {
            let mapping = aligner.map(s.as_slice(), true, false, None, None).unwrap_or_else(|err| panic!("Unable to align read with seq {s:?} to repeat-compressed reference for {repeat}\n{err}", s=s.to_ascii_uppercase()));
            for read in mapping {
                if let Some(s) = parse_cs(read, args.minlen, flanking, &repeat) {
                    // slice out inserted sequences from the CS tag
                    insertions.push(s.to_uppercase())
                }
            }
        }
        if insertions.len() < args.support {
            // Return a missing genotype if not enough insertions are found
            // this is too lenient - the support parameter is meant to be per haplotype
            return Ok(crate::vcf::VCFRecord::missing_genotype(
                &repeat,
                &repeat_ref_seq,
                insertions.len().to_string(),
            ));
        }
        debug!("{repeat}: Phasing {} insertions", insertions.len(),);
        let phased = crate::phase_insertions::split(&insertions, &repeat, args.find_outliers);
        match phased.hap2 {
            Some(phase2) => {
                consenses.push(crate::consensus::consensus(
                    &phased.hap1,
                    args.support,
                    &repeat,
                ));
                consenses.push(crate::consensus::consensus(&phase2, args.support, &repeat));
                if args.somatic {
                    // store all inserted sequences for identifying somatic variation
                    all_insertions.push(phased.hap1.join(":"));
                    all_insertions.push(phase2.join(":"));
                }
            }
            None => {
                // there was only one haplotype, homozygous, so this gets duplicated for reporting
                let consensus = crate::consensus::consensus(&phased.hap1, args.support, &repeat);
                consenses.push(consensus.clone());
                consenses.push(consensus);
                if args.somatic {
                    // store all inserted sequences for identifying somatic variation
                    all_insertions.push(phased.hap1.join(":"));
                }
                if let Some(splitflag) = phased.flag {
                    flags.push(splitflag);
                }
            }
        }
    } else {
        for phase in [1, 2] {
            let mut insertions = vec![];
            // get the sequences of this phase (or all if unphased)
            let seq = seqs.get(&phase).unwrap();
            debug!("{repeat}: Phase {}: Aliging {} reads", phase, seq.len());
            // align the reads to the new repeat-compressed reference
            for s in seq {
                let mapping = aligner.map(s.as_slice(), true, false, None, None).unwrap_or_else(|err| panic!("Unable to align read with seq {s:?} to repeat-compressed reference for {repeat}\n{err}", s=s.to_ascii_uppercase()));
                for read in mapping {
                    if let Some(s) = parse_cs(read, args.minlen, flanking, &repeat) {
                        // slice out inserted sequences from the CS tag
                        insertions.push(s.to_uppercase())
                    }
                }
            }
            debug!(
                "{repeat}: Phase {}: Creating consensus from {} insertions",
                phase,
                insertions.len(),
            );
            consenses.push(crate::consensus::consensus(
                &insertions,
                args.support,
                &repeat,
            ));

            if args.somatic {
                // store all inserted sequences for identifying somatic variation
                all_insertions.push(insertions.join(":"));
            }
        }
    }
    Ok(crate::vcf::VCFRecord::new(
        consenses,
        repeat_ref_seq,
        args.somatic,
        all_insertions,
        repeat,
        flags,
    ))
}

fn get_overlapping_reads(
    bamf: &String,
    repeat: &crate::repeats::RepeatInterval,
    unphased: bool,
) -> Option<HashMap<u8, Vec<Vec<u8>>>> {
    let mut bam = if bamf.starts_with("s3") || bamf.starts_with("https://") {
        bam::IndexedReader::from_url(&Url::parse(bamf.as_str()).expect("Failed to parse s3 URL"))
            .unwrap_or_else(|err| panic!("Error opening remote BAM: {err}"))
    } else {
        bam::IndexedReader::from_path(bamf)
            .unwrap_or_else(|err| panic!("Error opening local BAM: {err}"))
    };
    let tid = bam
        .header()
        .tid(repeat.chrom.as_bytes())
        .unwrap_or_else(|| panic!("Invalid chromosome {}", repeat.chrom));
    bam.fetch((tid, repeat.start, repeat.end))
        .unwrap_or_else(|err| panic!("Failure to extract reads from bam for {repeat}:\n{err}"));
    // Per haplotype the read sequences are kept in a dictionary
    let mut seqs = HashMap::from([(0, Vec::new()), (1, Vec::new()), (2, Vec::new())]);

    // extract sequences spanning the repeat locus
    for r in bam.rc_records() {
        let r = r.unwrap_or_else(|err| panic!("Error reading BAM file in region {repeat}:\n{err}"));
        if unphased {
            // if unphased put reads in phase 0
            seqs.get_mut(&0).unwrap().push(r.seq().as_bytes());
        } else {
            let phase = crate::utils::get_phase(&r);
            if phase > 0 {
                let seq = r.seq().as_bytes();
                seqs.get_mut(&phase).unwrap().push(seq);
                // writing fasta to stdout
                // println!(">read_{}\n{}", phase, std::str::from_utf8(&seq).unwrap());
            }
        }
    }
    if seqs.is_empty() {
        // error/warning message depends on whether we are looking for phased reads or not
        if unphased {
            eprintln!("Cannot genotype repeat at {repeat}: no reads found");
        } else {
            eprintln!("Cannot genotype repeat at {repeat}: no phased reads found");
        }
        None
    } else {
        Some(seqs)
    }
}

fn parse_cs(
    read: Mapping,
    minlen: usize,
    flanking: u32,
    repeat: &crate::repeats::RepeatInterval,
) -> Option<String> {
    // parses the CS tag of a <read> and returns the inserted sequence if it is longer than <minlen>
    // the reads are aligned to the repeat compressed reference genome,
    // which was constructed with <flanking> number of bases up and downstream of the repeat
    let mut ref_pos = read.target_start;
    let alignment = read.alignment.expect("Unable to access alignment");
    let cs = alignment.cs.expect("Unable to get the cs field");

    // split the cs tag using the regular expression
    // e.g. ':32*nt*na:10-gga:5+aaa:10' into (':32', '*nt', '*na', ':10', '-gga', ':5', '+aaa', ':10')
    let re = Regex::new(r"(:\d+)|(\*\w+)|(\+\w+)|(-\w+)").unwrap();

    let mut insertions = Vec::new();
    let interval_around_junction = flanking as i32 - 10..=flanking as i32 + 10;

    for cap in re.captures_iter(&cs) {
        let op = &cap[0].chars().next().unwrap();
        match op {
            ':' => {
                // match consumes the reference (and the read)
                ref_pos += cap[0][1..]
                    .parse::<i32>()
                    .expect("Unable to parse length from CS':' operation");
            }
            '*' => {
                // mismatch consumes the reference (and the read)
                // the cs tag is of the form *na, where n is the ref nucleotide and a the read nucleotide
                // therefore the number of nucleotides to advance is half the length of sequence from this cs operation
                ref_pos += cap[0][1..].len() as i32 / 2;
            }
            '-' => {
                // deletion consumes the reference
                // the cs tag is of the form -gga, where gga is the deleted sequence
                // therefore the number of nucleotides to advance is the length of sequence from this cs operation
                ref_pos += cap[0][1..].len() as i32;
            }
            '+' => {
                // insertion consumes the read
                // we only care about insertions that are within +/-10 bases of the repeat locus that was excised from the reference
                // this is the ref_pos
                // the cs tag is of the form +aaa, where aaa is the inserted sequence
                // if the insertion is longer than the minimum length, it is added to the list of insertions
                if cap[0][1..].len() > minlen && interval_around_junction.contains(&ref_pos) {
                    insertions.push(cap[0][1..].to_string());
                } else if cap[0][1..].len() > minlen {
                    debug!(
                        "{repeat}: Insertion {} is too far from the repeat locus junction to be considered: {}",
                        cap[0][1..].to_string(),
                        ref_pos
                    );
                }
            }
            _ => {
                // error
                panic!("Unexpected operation in cs tag: {op}");
            }
        }
    }
    if !insertions.is_empty() {
        Some(insertions.join(""))
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_get_overlapping_reads() {
        let bam = String::from("test_data/small-test-phased.bam");
        let repeat = crate::repeats::RepeatInterval {
            chrom: String::from("chr7"),
            start: 154654404,
            end: 154654432,
        };
        let unphased = false;
        let _reads = get_overlapping_reads(&bam, &repeat, unphased);
    }

    #[test]
    fn test_parse_cs() {
        let bam = String::from("test_data/small-test-phased.bam");
        let fasta = String::from("test_data/chr7.fa.gz");
        let repeat = crate::repeats::RepeatInterval {
            chrom: String::from("chr7"),
            start: 154654404,
            end: 154654432,
        };
        let flanking = 2000;
        let minlen = 5;
        let _support = 1;
        let unphased = false;
        let repeat_compressed_reference = repeat.make_repeat_compressed_sequence(&fasta, flanking);
        let binding = get_overlapping_reads(&bam, &repeat, unphased).unwrap();
        let read = binding
            .get(&1)
            .expect("Could not find phase 1")
            .iter()
            .next()
            .expect("No reads found");
        let aligner = minimap2::Aligner::builder()
            .map_ont()
            .with_cigar()
            .with_seq(&repeat_compressed_reference)
            .expect("Unable to build index");
        let mapping = aligner.map(read.as_slice(), true, false, None, None).unwrap_or_else(|_| panic!("Unable to align read with seq {read:?} to repeat-compressed reference for {repeat}", read=read.to_ascii_uppercase()));

        let _insertion = parse_cs(
            mapping
                .first()
                .expect("Could not get first mapping")
                .clone(),
            minlen,
            flanking,
            &repeat,
        );
    }

    #[test]
    fn test_genotype_repeat() {
        let repeat = crate::repeats::RepeatInterval {
            chrom: String::from("chr7"),
            start: 154654404,
            end: 154654432,
        };
        let args = Cli {
            bam: String::from("test_data/small-test-phased.bam"),
            fasta: String::from("test_data/chr7.fa.gz"),
            region: None,
            region_file: None,
            minlen: 5,
            support: 1,
            somatic: false,
            unphased: false,
            find_outliers: false,
            threads: 1,
            sample: None,
        };

        let genotype = genotype_repeat(repeat, &args);
        println!("{}", genotype.expect("Unable to genotype repeat"));
    }

    #[test]
    fn test_genotype_repeat_unphased() {
        let args = Cli {
            bam: String::from("test_data/small-test-phased.bam"),
            fasta: String::from("test_data/chr7.fa.gz"),
            region: None,
            region_file: None,
            minlen: 5,
            support: 1,
            somatic: false,
            unphased: true,
            find_outliers: false,
            threads: 1,
            sample: None,
        };
        let repeat = crate::repeats::RepeatInterval {
            chrom: String::from("chr7"),
            start: 154654404,
            end: 154654432,
        };

        let genotype = genotype_repeat(repeat, &args);
        println!("{}", genotype.expect("Unable to genotype repeat"));
    }

    #[test]
    fn test_genotype_repeat_somatic() {
        let args = Cli {
            bam: String::from("test_data/small-test-phased.bam"),
            fasta: String::from("test_data/chr7.fa.gz"),
            region: None,
            region_file: None,
            minlen: 5,
            support: 1,
            somatic: true,
            unphased: false,
            find_outliers: false,
            threads: 1,
            sample: None,
        };
        let repeat = crate::repeats::RepeatInterval {
            chrom: String::from("chr7"),
            start: 154654404,
            end: 154654432,
        };

        let genotype = genotype_repeat(repeat, &args);
        println!("{}", genotype.expect("Unable to genotype repeat"));
    }

    #[test]
    #[ignore]
    fn test_genotype_repeat_s3() {
        let args = Cli {
            // bam:  String::from("https://s3.amazonaws.com/1000g-ont/aligned_data_minimap2_2.24/HG01312/aligned_bams/HG01312.STD_eee-prom1_guppy-6.3.7-sup-prom_fastq_pass.phased.bam");
            // bam: String::from("s3://1000g-ont/aligned_data_minimap2_2.24/HG01312/aligned_bams/HG01312.STD_eee-prom1_guppy-6.3.7-sup-prom_fastq_pass.phased.bam");
            bam: String::from("https://1000g-ont.s3.amazonaws.com/aligned_data_minimap2_2.24/HG01312/aligned_bams/HG01312.STD_eee-prom1_guppy-6.3.7-sup-prom_fastq_pass.phased.bam"),
            fasta: String::from("test_data/chr7.fa.gz"),
            region: None,
            region_file: None,
            minlen: 5,
            support: 1,
            somatic: true,
            unphased: false,
            find_outliers: false,
            threads: 1,
            sample: None,
        };

        let repeat = crate::repeats::RepeatInterval {
            chrom: String::from("chr7"),
            start: 154654404,
            end: 154654432,
        };
        let genotype = genotype_repeat(repeat, &args);
        println!("{}", genotype.expect("Unable to genotype repeat"));
    }
}
