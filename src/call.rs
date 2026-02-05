use crate::repeats::RepeatIntervalIterator;
use indicatif::{ParallelProgressIterator, ProgressBar, ProgressIterator, ProgressStyle};
use log::{debug, error};
use rayon::prelude::*;
use std::cell::RefCell;
use std::io::Write;
use std::{io, sync::Mutex};

use crate::{Cli, genotype, parse_bam, utils::check_files_exist};

// Thread-local storage for BAM readers to avoid creating a new reader for each locus
// This dramatically reduces file descriptor usage when processing many loci in parallel
thread_local! {
    static BAM_READER: RefCell<Option<rust_htslib::bam::IndexedReader>> = const { RefCell::new(None) };
}

pub fn genotype_repeats(args: Cli) {
    debug!("Genotyping STRs in {}", args.bam);
    check_files_exist(&args);
    let repeats = get_targets(&args);
    crate::vcf::write_vcf_header(&args);
    let stdout = io::stdout(); // get the global stdout entity
    let mut handle = io::BufWriter::new(stdout); // wrap that handle in a buffer

    // Setup progress bar with smoothed ETA
    let num_intervals = repeats.num_intervals;
    let pb = ProgressBar::new(num_intervals as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} loci ({eta})")
            .expect("Failed to set progress bar template")
    );
    // Enable steady tick for smoothed ETA calculation (updates every 100ms)
    pb.enable_steady_tick(std::time::Duration::from_millis(100));

    if args.threads == 1 {
        // When running single threaded things become easier and the tool will require less memory
        // Output is returned in the same order as the bed, and therefore not sorted before writing immediately to stdout
        // The indexedreader is created once and passed on to the function
        let mut bam = parse_bam::create_bam_reader(&args.bam, &args.fasta);
        for mut repeat in repeats.progress_with(pb.clone()) {
            if let Ok(output) =
                genotype::genotype_repeat_singlethreaded(&mut repeat, &args, &mut bam)
            {
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
        // Pre-allocate with exact capacity to avoid reallocation
        let genotypes = Mutex::new(Vec::with_capacity(num_intervals));
        // par_bridge does not guarantee that results are returned in order
        repeats
            .par_bridge()
            .progress_with(pb.clone())
            .for_each(|mut repeat| {
                // Use thread-local BAM reader to avoid file descriptor exhaustion
                // Each thread maintains its own reader and reuses it for all loci
                BAM_READER.with(|reader_cell| {
                    let mut reader_opt = reader_cell.borrow_mut();
                    if reader_opt.is_none() {
                        // First use in this thread - create a new reader
                        *reader_opt = Some(parse_bam::create_bam_reader(&args.bam, &args.fasta));
                    }
                    let reader = reader_opt.as_mut().unwrap();
                    if let Ok(output) =
                        genotype::genotype_repeat_singlethreaded(&mut repeat, &args, reader)
                    {
                        let mut geno = genotypes.lock().expect("Unable to lock genotypes mutex");
                        geno.push(output);
                    } else {
                        error!("Problem processing {repeat}");
                    }
                });
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
        repeats
            .par_bridge()
            .progress_with(pb)
            .for_each(|mut repeat| {
                // Use thread-local BAM reader to avoid file descriptor exhaustion
                // Each thread maintains its own reader and reuses it for all loci
                BAM_READER.with(|reader_cell| {
                    let mut reader_opt = reader_cell.borrow_mut();
                    if reader_opt.is_none() {
                        // First use in this thread - create a new reader
                        *reader_opt = Some(parse_bam::create_bam_reader(&args.bam, &args.fasta));
                    }
                    let reader = reader_opt.as_mut().unwrap();
                    if let Ok(output) =
                        genotype::genotype_repeat_singlethreaded(&mut repeat, &args, reader)
                    {
                        println!("{output}");
                    }
                });
            });
    }
}

fn get_targets(args: &Cli) -> RepeatIntervalIterator {
    let repeats = match (&args.region, &args.region_file, args.pathogenic) {
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
    };

    // Filter by max_locus size if specified
    if let Some(max_size) = args.max_locus {
        repeats.filter_by_size(max_size)
    } else {
        repeats
    }
}
