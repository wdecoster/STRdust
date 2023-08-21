use bio::io::bed;
use log::{debug, error, info};
use rayon::prelude::*;
use regex::Regex;
use rust_htslib::{bam, bam::Read};
use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::Mutex;
use url::Url;

use minimap2::*;

pub fn genotype_repeats(
    bamp: PathBuf,
    fasta: PathBuf,
    region: Option<String>,
    region_file: Option<PathBuf>,
    minlen: usize,
    support: usize,
    threads: usize,
    sample: Option<String>,
    somatic: bool,
    unphased: bool,
) {
    if !bamp.is_file() {
        error!(
            "ERROR: path to bam file {} is not valid!\n\n",
            bamp.display()
        );
        panic!();
    };
    let bamf = bamp.into_os_string().into_string().unwrap();
    let fastaf = fasta.into_os_string().into_string().unwrap();
    debug!("Genotyping STRs in {bamf}");
    match (region, region_file) {
        (Some(_region), Some(_region_file)) => {
            error!("ERROR: Specify either a region (-r) or region_file (-R), not both!\n\n");
            panic!();
        }
        (None, None) => {
            error!("ERROR: Specify one of region (-r) or region_file (-R)!\n\n");
            panic!();
        }
        (Some(region), None) => {
            let (chrom, start, end) =
                crate::utils::process_region(region).expect("Error when parsing the region string");
            match genotype_repeat(
                &bamf, &fastaf, chrom, start, end, minlen, support, somatic, unphased,
            ) {
                Ok(output) => {
                    crate::utils::write_vcf_header(&fastaf, &bamf, sample);
                    println!("{output}")
                }
                Err(chrom) => error!("Contig {chrom} not found in bam file"),
            };
        }
        (None, Some(region_file)) => {
            // TODO: check if bed file is okay?
            let mut reader =
                bed::Reader::from_file(region_file.into_os_string().into_string().unwrap())
                    .expect("Problem reading bed file!");
            crate::utils::write_vcf_header(&fastaf, &bamf, sample);
            if threads > 1 {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .unwrap();
                // chrom_reported and genotypes are vectors that are used by multiple threads to add findings, therefore as a Mutex
                // chrom_reported contains those chromosomes for which an error (absence in the bam) was already reported
                // to avoid reporting the same error multiple times
                let chrom_reported = Mutex::new(Vec::new());
                // genotypes contains the output of the genotyping, a struct instance
                // let genotypes = Mutex::new(Vec::new());
                // par_bridge does not guarantee that results are returned in order
                reader.records().par_bridge().for_each(|record| {
                    let rec = record.expect("Error reading bed record.");
                    match genotype_repeat(
                        &bamf,
                        &fastaf,
                        rec.chrom().to_string(),
                        rec.start().try_into().unwrap(),
                        rec.end().try_into().unwrap(),
                        minlen,
                        support,
                        somatic,
                        unphased,
                    ) {
                        Ok(output) => {
                            // let mut geno = genotypes.lock().unwrap();
                            // geno.push(output);
                            println!("{output}")
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
                // let mut genotypes_vec = genotypes.lock().unwrap();
                // The final output is sorted by chrom, start and end
                // genotypes_vec.sort_unstable();
                // for g in &mut *genotypes_vec {
                //     println!("{g}");
                // }
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
                        &bamf,
                        &fastaf,
                        rec.chrom().to_string(),
                        rec.start().try_into().unwrap(),
                        rec.end().try_into().unwrap(),
                        minlen,
                        support,
                        somatic,
                        unphased,
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
    bamf: &String,
    fasta: &String,
    chrom: String,
    start: u32,
    end: u32,
    minlen: usize,
    support: usize,
    somatic: bool,
    unphased: bool,
) -> Result<String, String> {
    let flanking = 2000;
    let newref = crate::repeat_compressed_ref::make_repeat_compressed_sequence(
        fasta, &chrom, start, end, flanking,
    );
    if newref.is_none() {
        // Return a missing genotype if the repeat is not found in the fasta file
        return Ok(crate::write_vcf::missing_genotype(&chrom, start, end, "N"));
    }
    let (repeat_compressed_reference, repeat_ref_seq) = newref.unwrap();

    let reads = get_overlapping_reads(bamf, chrom.clone(), start, end, unphased);
    if reads.is_none() {
        // Return a missing genotype if no (phased) reads overlap the repeat
        return Ok(crate::write_vcf::missing_genotype(
            &chrom,
            start,
            end,
            &repeat_ref_seq,
        ));
    }
    let seqs = reads.unwrap();

    // Create an index for minimap2 alignment to the artificial reference
    let aligner = minimap2::Aligner::builder()
        .map_ont()
        .with_cigar()
        .with_seq(&repeat_compressed_reference)
        .unwrap_or_else(|err| panic!("Unable to build index:\n{err}"));
    let mut consenses: Vec<Result<(String, usize, usize), crate::consensus::ConsensusError>> =
        vec![];
    let mut all_insertions = vec![]; // only used with `--somatic`
    if unphased {
        let mut insertions = vec![];
        // get the sequences of this phase (or all if unphased)
        let seq = seqs.get(&0).unwrap();
        debug!("Unphased: Aliging {} reads", seq.len());
        // align the reads to the new repeat-compressed reference
        for s in seq {
            let mapping = aligner.map(s.as_slice(), true, false, None, None).unwrap_or_else(|err| panic!("Unable to align read with seq {s:?} to repeat-compressed reference for {chrom}:{start}-{end}\n{err}", s=s.to_ascii_uppercase(), chrom=chrom, start=start, end=end));
            for read in mapping {
                if let Some(s) = parse_cs(read, minlen, flanking) {
                    // slice out inserted sequences from the CS tag
                    insertions.push(s.to_uppercase())
                }
            }
        }
        if insertions.len() < support {
            // Return a missing genotype if not enough insertions are found
            return Ok(crate::write_vcf::missing_genotype(
                &chrom,
                start,
                end,
                &repeat_ref_seq,
            ));
        }
        debug!("Splitting {} insertions in two phases", insertions.len(),);
        let phased = crate::phase_insertions::split(&insertions);
        match phased {
            (Some(phase1), Some(phase2)) => {
                consenses.push(crate::consensus::consensus(&phase1, support));
                consenses.push(crate::consensus::consensus(&phase2, support));
                if somatic {
                    // store all inserted sequences for identifying somatic variation
                    all_insertions.push(phase1.join(","));
                    all_insertions.push(phase2.join(","));
                }
            }
            (Some(phase1), None) => {
                let consensus = crate::consensus::consensus(&phase1, support);
                consenses.push(consensus.clone());
                consenses.push(consensus);
                if somatic {
                    // store all inserted sequences for identifying somatic variation
                    all_insertions.push(phase1.join(","));
                }
            }
            _ => {
                error!("Unexpected scenario after haplotype splitting");
                panic!();
            }
        }
    } else {
        for phase in [1, 2] {
            let mut insertions = vec![];
            // get the sequences of this phase (or all if unphased)
            let seq = seqs.get(&phase).unwrap();
            debug!("Phase {}: Aliging {} reads", phase, seq.len());
            // align the reads to the new repeat-compressed reference
            for s in seq {
                let mapping = aligner.map(s.as_slice(), true, false, None, None).unwrap_or_else(|err| panic!("Unable to align read with seq {s:?} to repeat-compressed reference for {chrom}:{start}-{end}\n{err}", s=s.to_ascii_uppercase(), chrom=chrom, start=start, end=end));
                for read in mapping {
                    if let Some(s) = parse_cs(read, minlen, flanking) {
                        // slice out inserted sequences from the CS tag
                        insertions.push(s.to_uppercase())
                    }
                }
            }
            debug!(
                "Phase {}: Creating consensus from {} insertions",
                phase,
                insertions.len(),
            );
            consenses.push(crate::consensus::consensus(&insertions, support));

            if somatic {
                // store all inserted sequences for identifying somatic variation
                all_insertions.push(insertions.join(","));
            }
        }
    }
    Ok(crate::write_vcf::write_vcf(
        consenses,
        repeat_ref_seq,
        somatic,
        all_insertions,
        chrom,
        start,
        end,
    ))
}

fn get_overlapping_reads(
    bamf: &String,
    chrom: String,
    start: u32,
    end: u32,
    unphased: bool,
) -> Option<HashMap<u8, Vec<Vec<u8>>>> {
    let mut bam = if bamf.starts_with("s3") || bamf.starts_with("https://") {
        bam::IndexedReader::from_url(&Url::parse(bamf).expect("Failed to parse s3 URL"))
            .unwrap_or_else(|err| panic!("Error opening remote BAM: {err}"))
    } else {
        bam::IndexedReader::from_path(bamf)
            .unwrap_or_else(|err| panic!("Error opening local BAM: {err}"))
    };
    info!("Checks passed, genotyping repeat");
    let tid = bam
        .header()
        .tid(chrom.as_bytes())
        .unwrap_or_else(|| panic!("Invalid chromosome {}", chrom));
    bam.fetch((tid, start, end)).unwrap_or_else(|err| {
        panic!("Failure to extract reads from bam for {chrom}:{start}-{end}:\n{err}")
    });
    // Per haplotype the read sequences are kept in a dictionary
    let mut seqs = HashMap::from([(0, Vec::new()), (1, Vec::new()), (2, Vec::new())]);

    // extract sequences spanning the repeat locus
    for r in bam.rc_records() {
        let r = r.unwrap_or_else(|err| {
            panic!("Error reading BAM file in region {chrom}:{start}-{end}:\n{err}")
        });
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
            eprintln!("Cannot genotype repeat at {chrom}:{start}-{end}: no reads found");
        } else {
            eprintln!("Cannot genotype repeat at {chrom}:{start}-{end}: no phased reads found");
        }
        None
    } else {
        Some(seqs)
    }
}

fn parse_cs(read: Mapping, minlen: usize, flanking: u32) -> Option<String> {
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
                if cap[0][1..].len() > minlen
                    && (flanking as i32 - 10..=flanking as i32 + 10).contains(&ref_pos)
                {
                    insertions.push(cap[0][1..].to_string());
                } else if cap[0][1..].len() > minlen {
                    debug!(
                        "Insertion {} is too far from the repeat locus junction to be considered: {}",
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
        let chrom = String::from("chr7");
        let start = 154654404;
        let end = 154654432;
        let unphased = false;
        let _reads = get_overlapping_reads(&bam, chrom, start, end, unphased);
    }

    #[test]
    fn test_parse_cs() {
        let bam = String::from("test_data/small-test-phased.bam");
        let fasta = String::from("test_data/chr7.fa.gz");
        let chrom = String::from("chr7");
        let start = 154654404;
        let end = 154654432;
        let flanking = 2000;
        let minlen = 5;
        let _support = 1;
        let unphased = false;
        let (repeat_compressed_reference, _) =
            crate::repeat_compressed_ref::make_repeat_compressed_sequence(
                &fasta, &chrom, start, end, flanking,
            )
            .expect("Unable to make repeat compressed sequence");
        let binding = get_overlapping_reads(&bam, chrom.clone(), start, end, unphased).unwrap();
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
        let mapping = aligner.map(read.as_slice(), true, false, None, None).unwrap_or_else(|_| panic!("Unable to align read with seq {read:?} to repeat-compressed reference for {chrom}:{start}-{end}", read=read.to_ascii_uppercase(), chrom=chrom, start=start, end=end));

        let _insertion = parse_cs(
            mapping
                .first()
                .expect("Could not get first mapping")
                .clone(),
            minlen,
            flanking,
        );
    }

    #[test]
    fn test_genotype_repeat() {
        let bam = String::from("test_data/small-test-phased.bam");
        let fasta = String::from("test_data/chr7.fa.gz");
        let chrom = String::from("chr7");
        let start = 154654404;
        let end = 154654432;
        let minlen = 5;
        let support = 1;
        let somatic = false;
        let unphased = false;
        let genotype = genotype_repeat(
            &bam, &fasta, chrom, start, end, minlen, support, somatic, unphased,
        );
        println!("{}", genotype.expect("Unable to genotype repeat"));
    }

    #[test]
    fn test_genotype_repeat_unphased() {
        let bam = String::from("test_data/small-test-phased.bam");
        let fasta = String::from("test_data/chr7.fa.gz");
        let chrom = String::from("chr7");
        let start = 154654404;
        let end = 154654432;
        let minlen = 5;
        let support = 1;
        let somatic = false;
        let unphased = true;
        let genotype = genotype_repeat(
            &bam, &fasta, chrom, start, end, minlen, support, somatic, unphased,
        );
        println!("{}", genotype.expect("Unable to genotype repeat"));
    }

    #[test]
    fn test_genotype_repeat_somatic() {
        let bam = String::from("test_data/small-test-phased.bam");
        let fasta = String::from("test_data/chr7.fa.gz");
        let chrom = String::from("chr7");
        let start = 154654404;
        let end = 154654432;
        let minlen = 5;
        let support = 1;
        let somatic = true;
        let unphased = false;
        let genotype = genotype_repeat(
            &bam, &fasta, chrom, start, end, minlen, support, somatic, unphased,
        );
        println!("{}", genotype.expect("Unable to genotype repeat"));
    }

    #[test]
    #[ignore]
    fn test_genotype_repeat_s3() {
        // let bam = String::from("https://s3.amazonaws.com/1000g-ont/aligned_data_minimap2_2.24/HG01312/aligned_bams/HG01312.STD_eee-prom1_guppy-6.3.7-sup-prom_fastq_pass.phased.bam");
        let bam = String::from("s3://1000g-ont/aligned_data_minimap2_2.24/HG01312/aligned_bams/HG01312.STD_eee-prom1_guppy-6.3.7-sup-prom_fastq_pass.phased.bam");
        let fasta = String::from("test_data/chr7.fa.gz");
        let chrom = String::from("chr7");
        let start = 154654404;
        let end = 154654432;
        let minlen = 5;
        let support = 1;
        let somatic = true;
        let unphased = false;
        let genotype = genotype_repeat(
            &bam, &fasta, chrom, start, end, minlen, support, somatic, unphased,
        );
        println!("{}", genotype.expect("Unable to genotype repeat"));
    }
}
