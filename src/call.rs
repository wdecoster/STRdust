use bio::io::bed;
use log::{debug, error};
use rayon::prelude::*;
use std::io::Write;
use std::{io, sync::Mutex};

use crate::{genotype, parse_bam, Cli};

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
                if let Ok(output) =
                    genotype::genotype_repeat_singlethreaded(&repeat, &args, &mut bam)
                {
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
            let stdout = io::stdout(); // get the global stdout entity
            let mut handle = io::BufWriter::new(stdout); // optional: wrap that handle in a buffer

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
                        if let Ok(output) = genotype::genotype_repeat_multithreaded(&repeat, &args)
                        {
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
                    writeln!(handle, "{g}").expect("Failed writing the result.");
                }
            } else {
                // When running single threaded things become easier and the tool will require less memory
                // Output is returned in the same order as the bed, and therefore not sorted before writing immediately to stdout
                // The indexedreader is created once and passed on to the function
                let mut bam = parse_bam::create_bam_reader(&args.bam);
                // genotypes contains the output of the genotyping, a struct instance
                for record in reader.records() {
                    let rec = record.expect("Error reading bed record.");
                    let repeat = crate::repeats::RepeatInterval::from_bed(&rec, &args.fasta);
                    if let Some(repeat) = repeat {
                        if let Ok(output) =
                            genotype::genotype_repeat_singlethreaded(&repeat, &args, &mut bam)
                        {
                            writeln!(handle, "{output}").expect("Failed writing the result.");
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
// Alternatively, this error can be raised if the end is smaller than the start
fn error_invalid_interval(rec: &bed::Record) {
    error!(
        "ERROR: Invalid interval in bed file: {chrom}:{begin}-{end}",
        chrom = rec.chrom(),
        begin = rec.start(),
        end = rec.end()
    );
}
