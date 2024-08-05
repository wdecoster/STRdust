use crate::repeats::RepeatIntervalIterator;
use indicatif::ParallelProgressIterator;
use indicatif::ProgressIterator;
use log::{debug, error};
use rayon::prelude::*;
use std::io::Write;
use std::{io, sync::Mutex};

use crate::{genotype, parse_bam, Cli, utils::check_files_exist};

pub fn genotype_repeats(args: Cli) {
    debug!("Genotyping STRs in {}", args.bam);
    check_files_exist(&args);
    let repeats = get_targets(&args);
    crate::vcf::write_vcf_header(&args);
    let stdout = io::stdout(); // get the global stdout entity
    let mut handle = io::BufWriter::new(stdout); // wrap that handle in a buffer
    if args.threads == 1 {
        // When running single threaded things become easier and the tool will require less memory
        // Output is returned in the same order as the bed, and therefore not sorted before writing immediately to stdout
        // The indexedreader is created once and passed on to the function
        let num_intervals = repeats.num_intervals;
        let mut bam = parse_bam::create_bam_reader(&args.bam, &args.fasta);
        for mut repeat in repeats.progress_count(num_intervals as u64) {
            if let Ok(output) = genotype::genotype_repeat_singlethreaded(&mut repeat, &args, &mut bam) {
                writeln!(handle, "{output}").expect("Failed writing the result.");
            }
        }
    } else if args.sorted {
        // output is sorted by chrom, start and end
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .expect("Failed to create threadpool");
        // genotypes contains the output of the genotyping, a struct instance
        let genotypes = Mutex::new(Vec::new());
        // par_bridge does not guarantee that results are returned in order
        let num_intervals = repeats.num_intervals;
        repeats
            .par_bridge()
            .progress_count(num_intervals as u64)
            .for_each(|mut repeat| {
                if let Ok(output) = genotype::genotype_repeat_multithreaded(&mut repeat, &args) {
                    let mut geno = genotypes.lock().expect("Unable to lock genotypes mutex");
                    geno.push(output);
                } else {
                    error!("Problem processing {repeat}");
                }
            });
        let mut genotypes_vec = genotypes.lock().unwrap();
        // The final output is sorted by chrom, start and end
        genotypes_vec.sort_unstable();
        for g in &mut *genotypes_vec {
            writeln!(handle, "{g}").expect("Failed writing the result.");
        }
    } else {
        // run in parallel, but output is most probably unsorted
        // par_bridge does not guarantee that results are returned in order
        // but writing it out immediately makes for low memory usage
        // this does not use the BufWriter (as that doesn't work across threads), but writes directly to stdout
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .expect("Failed to create threadpool");
        let num_intervals = repeats.num_intervals;
        repeats
            .par_bridge()
            .progress_count(num_intervals as u64)
            .for_each(|mut repeat| {
                if let Ok(output) = genotype::genotype_repeat_multithreaded(&mut repeat, &args) {
                    println!("{output}");
                }
            });
    }
}

fn get_targets(args: &Cli) -> RepeatIntervalIterator {
    match (&args.region, &args.region_file, args.pathogenic) {
        // a region string
        (Some(region), None, false) => RepeatIntervalIterator::from_string(region, &args.fasta),
        // a region file
        (None, Some(region_file), false) => {
            RepeatIntervalIterator::from_bed(region_file, &args.fasta)
        }
        // with --pathogenic
        (None, None, true) => RepeatIntervalIterator::pathogenic(&args.fasta),
        // invalid input
        _ => {
            eprintln!("ERROR: Specify a region string (-r), a region_file (-R) or --pathogenic!\n");
            std::process::exit(1);
        }
    }
}
