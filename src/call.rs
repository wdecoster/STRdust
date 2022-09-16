use bio::io::bed;
use human_sort::compare as human_compare;
use log::{error, info};
use rayon::prelude::*;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::{Aux, Cigar};
use rust_htslib::faidx;
use rust_htslib::{bam, bam::Read};
use std::cmp::Ordering;
use std::collections::HashMap;
use std::f64::NAN;
use std::fmt;
use std::path::PathBuf;
use std::rc::Rc;
use std::sync::Mutex;

pub fn genotype_repeats(
    bamp: PathBuf,
    fasta: PathBuf,
    region: Option<String>,
    region_file: Option<PathBuf>,
    minlen: u32,
    support: usize,
    threads: usize,
) {
    if !bamp.is_file() {
        error!(
            "ERROR: path to bam file {} is not valid!\n\n",
            bamp.display()
        );
        panic!();
    };
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
            let (chrom, start, end) = crate::utils::process_region(region).unwrap();
            let bamf = bamp.into_os_string().into_string().unwrap();
            let fastaf = fasta.into_os_string().into_string().unwrap();
            match genotype_repeat(&bamf, &fastaf, chrom, start, end, minlen, support) {
                Ok(output) => println!("{}", output),
                Err(chrom) => error!("Contig {chrom} not found in bam file"),
            };
        }
        (None, Some(region_file)) => {
            if threads > 1 {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .unwrap();
                // TODO: check if bed file is okay
                let mut reader =
                    bed::Reader::from_file(region_file.into_os_string().into_string().unwrap())
                        .unwrap();
                let bamf = bamp.into_os_string().into_string().unwrap();
                let fastaf = fasta.into_os_string().into_string().unwrap();
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
                        &bamf,
                        &fastaf,
                        rec.chrom().to_string(),
                        rec.start().try_into().unwrap(),
                        rec.end().try_into().unwrap(),
                        minlen,
                        support,
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
                    println!("{}", g);
                }
            } else {
                // When running single threaded things become easier and the tool will require less memory
                // Output is returned in the same order as the bed, and therefore not sorted before writing to stdout
                // TODO: check if bed file is okay
                let mut reader =
                    bed::Reader::from_file(region_file.into_os_string().into_string().unwrap())
                        .unwrap();
                let bamf = bamp.into_os_string().into_string().unwrap();
                let fastaf = fasta.into_os_string().into_string().unwrap();
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
                    ) {
                        Ok(output) => {
                            println!("{}", output);
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
    _minlen: u32,
    _support: usize,
) -> Result<String, String> {
    let newref = make_repeat_compressed_sequence(fasta, &chrom, start, end);
    let mut bam = match bam::IndexedReader::from_path(&bamf) {
        Ok(handle) => handle,
        Err(e) => {
            error!("Error opening BAM {}.\n{}", bamf, e);
            panic!();
        }
    };
    info!("Checks passed, genotyping repeat");
    if let Some(tid) = bam.header().tid(chrom.as_bytes()) {
        bam.fetch((tid, start, end)).unwrap();
        // Per haplotype the read sequences are kept in a dictionary
        let mut seqs: HashMap<u8, Vec<Vec<u8>>> = HashMap::from([(1, Vec::new()), (2, Vec::new())]);

        // extract sequences spanning the repeat
        for r in bam.rc_records() {
            let r = r.expect("Error reading BAM file in region {chrom}:{start}-{end}.");
            let phase = get_phase(&r);
            if phase > 0 {
                let seq = r.seq().as_bytes();
                seqs.get_mut(&phase).unwrap().push(seq);
            }
        }
    }
    Ok("hiyaaa".to_string())
}

fn make_repeat_compressed_sequence(fasta: &String, chrom: &String, start: u32, end: u32) -> String {
    let fas = faidx::Reader::from_path(fasta).expect("Failed to read fasta");
    // TODO make sure we cannot try to read past the end of chromosomes
    let fas_left = fas
        .fetch_seq(&chrom, (start - 10000) as usize, start as usize)
        .expect("Failed to extract sequence from fasta");
    let fas_right = fas
        .fetch_seq(&chrom, end as usize, (end + 10000) as usize)
        .expect("Failed to extract sequence from fasta");
    let mut seq = std::str::from_utf8(fas_left)
        .expect("Failed to extract sequence from fasta")
        .to_string();
    seq.push_str(std::str::from_utf8(fas_right).expect("Failed to extract sequence from fasta"));
    seq
}

fn get_phase(record: &bam::Record) -> u8 {
    match record.aux(b"HP") {
        Ok(value) => {
            if let Aux::U8(v) = value {
                v
            } else {
                panic!("Unexpected type of Aux {:?}", value)
            }
        }
        Err(_e) => 0,
    }
}
