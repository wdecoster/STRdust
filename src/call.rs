use bio::io::bed;
use log::{debug, error};
use minimap2::*;
use rayon::prelude::*;
use regex::Regex;
use rust_htslib::bam;
use std::sync::Mutex;

use crate::{parse_bam, Cli};

pub fn genotype_repeats(args: Cli) {
    debug!("Genotyping STRs in {}", args.bam);
    match (&args.region, &args.region_file) {
        // both input options are specified: invalid
        (Some(_region), Some(_region_file)) => {
            error!("ERROR: Specify either a region (-r) or region_file (-R), not both!\n\n");
            panic!();
        }
        // no input options are specified: invalid
        (None, None) => {
            error!("ERROR: Specify one of region (-r) or region_file (-R)!\n\n");
            panic!();
        }
        // a region string is specified: single threaded
        (Some(region), None) => {
            let repeat = crate::repeats::RepeatInterval::from_string(region, &args.fasta);
            let mut bam = parse_bam::create_bam_reader(&args.bam);
            if let Some(repeat) = repeat {
                if let Ok(output) = genotype_repeat_singlethreaded(&repeat, &args, &mut bam) {
                    println!("{output}");
                }
            } else {
                // None is returned when the interval from the bed file does not appear in the fasta file
                error!("Interval {region} not found in the fasta file, ignoring...");
            }
        }
        // a region file is specified: single or multithreaded
        (None, Some(region_file)) => {
            // TODO: check if bed file is okay?
            let mut reader =
                bed::Reader::from_file(region_file).expect("Problem reading bed file!");
            crate::vcf::write_vcf_header(&args.fasta, &args.bam, &args.sample);
            // if multithreaded, create a threadpool
            if args.threads > 1 {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(args.threads)
                    .build()
                    .expect("Failed to create threadpool");
                // genotypes contains the output of the genotyping, a struct instance
                let genotypes = Mutex::new(Vec::new());
                // par_bridge does not guarantee that results are returned in order
                reader.records().par_bridge().for_each(|record| {
                    let rec: bed::Record = record.expect("Error reading bed record.");
                    let repeat = crate::repeats::RepeatInterval::from_bed(&rec, &args.fasta);
                    if let Some(repeat) = repeat {
                        if let Ok(output) = genotype_repeat_multithreaded(&repeat, &args) {
                            let mut geno =
                                genotypes.lock().expect("Unable to lock genotypes mutex");
                            geno.push(output);
                        } else {
                            error!("Problem processing {repeat}");
                        }
                    } else {
                        error_invalid_interval(&rec);
                    }
                });
                let mut genotypes_vec = genotypes.lock().unwrap();
                // The final output is sorted by chrom, start and end
                genotypes_vec.sort_unstable();
                for g in &mut *genotypes_vec {
                    println!("{g}");
                }
            } else {
                // When running single threaded things become easier and the tool will require less memory
                // Output is returned in the same order as the bed, and therefore not sorted before writing immediately to stdout
                let mut bam = parse_bam::create_bam_reader(&args.bam);
                // genotypes contains the output of the genotyping, a struct instance
                for record in reader.records() {
                    let rec = record.expect("Error reading bed record.");
                    let repeat = crate::repeats::RepeatInterval::from_bed(&rec, &args.fasta);
                    if let Some(repeat) = repeat {
                        if let Ok(output) = genotype_repeat_singlethreaded(&repeat, &args, &mut bam)
                        {
                            println!("{output}");
                        }
                    } else {
                        error_invalid_interval(&rec);
                    }
                }
            }
        }
    }
}
// None is returned when the interval from the bed file does not appear in the fasta file
fn error_invalid_interval(rec: &bed::Record) {
    error!(
        "ERROR: Invalid interval in bed file: {chrom}:{begin}-{end}",
        chrom = rec.chrom(),
        begin = rec.start(),
        end = rec.end()
    );
}

// when running multithreaded, the indexedreader has to be created every time again
fn genotype_repeat_multithreaded(
    repeat: &crate::repeats::RepeatInterval,
    args: &Cli,
) -> Result<crate::vcf::VCFRecord, String> {
    let mut bam = parse_bam::create_bam_reader(&args.bam);
    genotype_repeat(repeat, args, &mut bam)
}

// when running singlethreaded, the indexedreader is created once and simply passed on
fn genotype_repeat_singlethreaded(
    repeat: &crate::repeats::RepeatInterval,
    args: &Cli,
    bam: &mut bam::IndexedReader,
) -> Result<crate::vcf::VCFRecord, String> {
    genotype_repeat(repeat, args, bam)
}

/// This function genotypes a particular repeat defined by chrom, start and end in the specified bam file
/// All indel cigar operations longer than minlen are considered
/// The bam file is expected to be phased using the HP tag
fn genotype_repeat(
    repeat: &crate::repeats::RepeatInterval,
    args: &Cli,
    bam: &mut bam::IndexedReader,
) -> Result<crate::vcf::VCFRecord, String> {
    let flanking = 2000;
    let mut flags = vec![];
    let repeat_ref_seq = match repeat.reference_repeat_sequence(&args.fasta) {
        Some(seq) => seq,
        // Return a missing genotype if the repeat is not found in the fasta file
        None => {
            return Ok(crate::vcf::VCFRecord::missing_genotype(
                repeat,
                "N",
                ".".to_string(),
            ))
        }
    };

    let repeat_compressed_reference = repeat.make_repeat_compressed_sequence(&args.fasta, flanking);

    let reads = match crate::parse_bam::get_overlapping_reads(bam, repeat, args.unphased) {
        Some(seqs) => seqs,
        None => {
            // Return a missing genotype if no (phased) reads overlap the repeat
            return Ok(crate::vcf::VCFRecord::missing_genotype(
                repeat,
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
    // only used with `--somatic`
    let mut all_insertions = if args.somatic { Some(vec![]) } else { None };
    // only used with `--find_outliers`
    let mut outlier_insertions = if args.find_outliers {
        Some(vec![])
    } else {
        None
    };
    if args.unphased {
        // get the sequences
        let seq = reads.seqs.get(&0).unwrap();
        debug!("{repeat}: Unphased: Aliging {} reads", seq.len());
        // align the reads to the new repeat-compressed reference
        let insertions = find_insertions(seq, &aligner, args.minlen, flanking, repeat);
        if insertions.len() < args.support {
            // Return a missing genotype if not enough insertions are found
            // this is too lenient - the support parameter is meant to be per haplotype
            return Ok(crate::vcf::VCFRecord::missing_genotype(
                repeat,
                &repeat_ref_seq,
                insertions.len().to_string(),
            ));
        }
        debug!("{repeat}: Phasing {} insertions", insertions.len(),);
        let phased = crate::phase_insertions::split(&insertions, repeat, args.find_outliers);
        match phased.hap2 {
            Some(phase2) => {
                consenses.push(crate::consensus::consensus(
                    &phased.hap1,
                    args.support,
                    repeat,
                ));
                consenses.push(crate::consensus::consensus(&phase2, args.support, repeat));
                // store all inserted sequences for identifying somatic variation
                if let Some(ref mut all_ins) = all_insertions {
                    all_ins.extend([phased.hap1.join(":"), phase2.join(":")]);
                }
            }
            None => {
                // there was only one haplotype, homozygous, so this gets duplicated for reporting
                // not sure if cloning is the best approach here, but this is only the case for unphased data
                // and therefore is typically for small datasets obtained through capture methods
                let consensus = crate::consensus::consensus(&phased.hap1, args.support, repeat);
                consenses.push(consensus.clone());
                consenses.push(consensus);
                // store all inserted sequences for identifying somatic variation
                if let Some(ref mut all_ins) = all_insertions {
                    all_ins.push(phased.hap1.join(":"));
                }
                // if looking for outliers, and outliers were found, store them

                // escalate the flag to the VCF
                if let Some(splitflag) = phased.flag {
                    flags.push(splitflag);
                }
            }
        }
        if let Some(ref mut outliers_vec) = outlier_insertions {
            if let Some(outliers_found) = phased.outliers {
                outliers_vec.push(outliers_found.join(","));
            }
        }
    } else {
        // input reads are phased
        for phase in [1, 2] {
            // get the sequences of this phase
            let seq = reads.seqs.get(&phase).unwrap();
            debug!("{repeat}: Phase {}: Aliging {} reads", phase, seq.len());
            let insertions = find_insertions(seq, &aligner, args.minlen, flanking, repeat);

            debug!(
                "{repeat}: Phase {}: Creating consensus from {} insertions",
                phase,
                insertions.len(),
            );
            consenses.push(crate::consensus::consensus(
                &insertions,
                args.support,
                repeat,
            ));

            if let Some(ref mut all_ins) = all_insertions {
                // store all inserted sequences for identifying somatic variation
                all_ins.push(insertions.join(":"));
            }
        }
    }
    Ok(crate::vcf::VCFRecord::new(
        consenses,
        repeat_ref_seq,
        all_insertions,
        outlier_insertions,
        repeat,
        reads.ps,
        flags,
    ))
}

// may adapt the function below to allow for multiple alignment methods later
fn find_insertions(
    seq: &Vec<Vec<u8>>,
    aligner: &Aligner,
    minlen: usize,
    flanking: u32,
    repeat: &crate::repeats::RepeatInterval,
) -> Vec<String> {
    let mut insertions = vec![];

    // align the reads to the new repeat-compressed reference
    for s in seq {
        let mapping = aligner.map(s.as_slice(), true, false, None, None).unwrap_or_else(|err| panic!("Unable to align read with seq {s:?} to repeat-compressed reference for {repeat}\n{err}", s=s.to_ascii_uppercase()));
        for read in mapping {
            if let Some(s) = parse_cs(read, minlen, flanking, repeat) {
                // slice out inserted sequences from the CS tag
                insertions.push(s.to_uppercase())
            }
        }
    }
    insertions
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
    let interval_around_junction = flanking as i32 - 15..=flanking as i32 + 15;

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
                        "{repeat}: Insertion at {} is too far from the junction to be considered: {}",
                        ref_pos,
                        cap[0][1..].to_string(),
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
        let mut bam = parse_bam::create_bam_reader(&bam);
        let binding = crate::parse_bam::get_overlapping_reads(&mut bam, &repeat, unphased).unwrap();
        let read = binding
            .seqs
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
        let mut bam = parse_bam::create_bam_reader(&args.bam);
        let genotype = genotype_repeat(&repeat, &args, &mut bam);
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
        let mut bam = parse_bam::create_bam_reader(&args.bam);
        let genotype = genotype_repeat(&repeat, &args, &mut bam);
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
        let mut bam = parse_bam::create_bam_reader(&args.bam);
        let genotype = genotype_repeat(&repeat, &args, &mut bam);
        println!("{}", genotype.expect("Unable to genotype repeat"));
    }

    #[test]
    fn test_genotype_repeat_url() {
        let args = Cli {
            bam: String::from("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/hg38/HG00096.hg38.cram"),
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
        let mut bam = parse_bam::create_bam_reader(&args.bam);
        let genotype = genotype_repeat(&repeat, &args, &mut bam);
        println!("{}", genotype.expect("Unable to genotype repeat"));
    }
}
