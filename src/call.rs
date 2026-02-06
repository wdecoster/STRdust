use crate::bam_pool::BamReaderPool;
use crate::batching::{Batch, create_batches};
use crate::repeats::RepeatIntervalIterator;
use indicatif::{ProgressBar, ProgressStyle};
use log::{debug, error, info};
use rayon::prelude::*;
use rust_htslib::bam::{self, Read, ext::BamRecordExtensions};
use std::io::Write;
use std::{io, sync::Mutex};

use crate::{Cli, genotype, parse_bam, utils::check_files_exist};

/// Process a batch of nearby STR targets with optimized single-fetch approach
///
/// This function:
/// 1. Fetches all reads in the batch region once
/// 2. Filters to only store reads overlapping at least one target
/// 3. Performs QUICKREF screening (CIGAR-only check) for reference calls
/// 4. Extracts sequences only for targets that need full genotyping
///
/// Returns vector of VCF records in the same order as batch.repeats
fn process_batch(
    batch: &Batch,
    bam: &mut bam::IndexedReader,
    args: &Cli,
) -> Result<Vec<crate::vcf::VCFRecord>, String> {
    // Get chromosome tid and do fetch in one scope
    let tid = match bam.header().tid(batch.chromosome.as_bytes()) {
        Some(t) => t,
        None => return Err(format!("Invalid chromosome {}", batch.chromosome)),
    };

    // Single fetch for entire batch region
    bam.fetch((tid, batch.start, batch.end)).map_err(|err| {
        format!(
            "Failure to extract reads from bam for batch {}:{}-{}:\n{err}",
            batch.chromosome, batch.start, batch.end
        )
    })?;

    // Track target information and only store overlapping records
    struct TargetInfo {
        has_variation: bool,
        record_indices: Vec<usize>,
    }

    let mut filtered_records = Vec::new();
    let mut target_info: std::collections::HashMap<usize, TargetInfo> =
        std::collections::HashMap::new();

    // Single pass: check overlaps, do CIGAR check, store only relevant records
    for record_result in bam.rc_records() {
        let record_rc = record_result.map_err(|err| {
            format!(
                "Error reading BAM file in batch {}:{}-{}:\n{err}",
                batch.chromosome, batch.start, batch.end
            )
        })?;

        let read_start = record_rc.reference_start() as u32;
        let read_end = record_rc.reference_end() as u32;

        let mut overlaps_any = false;
        let future_idx = filtered_records.len();

        // Check all targets for overlap
        // Repeats are sorted by start position, so we can break early
        for (target_idx, target) in batch.repeats.iter().enumerate() {
            // Early exit: if repeat starts after read ends, no subsequent repeats will overlap
            if target.start >= read_end {
                break;
            }

            if read_start < target.end && read_end > target.start {
                overlaps_any = true;

                let info = target_info.entry(target_idx).or_insert_with(|| TargetInfo {
                    has_variation: false,
                    record_indices: Vec::new(),
                });

                info.record_indices.push(future_idx);

                // QUICKREF: CIGAR check (only if not already found variation and alignment_all is not set)
                // When --alignment-all is set, skip QUICKREF optimization to force full alignment
                if !args.alignment_all && !info.has_variation {
                    let diff = parse_bam::calculate_all_length_diff_from_cigar(
                        &record_rc,
                        target.start,
                        target.end,
                    );
                    if diff != 0 {
                        info.has_variation = true;
                    }
                }
            }
        }

        // Only store record if it overlaps at least one target
        // Convert Rc<Record> to owned Record for storage
        if overlaps_any {
            filtered_records.push((*record_rc).clone());
        }
    }

    // Process each target based on what we learned
    let mut results = Vec::new();

    for (target_idx, repeat) in batch.repeats.iter().enumerate() {
        match target_info.get(&target_idx) {
            None => {
                // No overlapping reads - use existing function to handle this case
                let mut repeat_mut = repeat.clone();
                if args.debug {
                    repeat_mut.set_time_stamp();
                }
                match genotype::genotype_repeat_singlethreaded(&mut repeat_mut, args, bam) {
                    Ok(vcf_record) => results.push(vcf_record),
                    Err(e) => {
                        error!("Problem processing {repeat}: {e}");
                    }
                }
            }
            Some(info) if !args.alignment_all && !info.has_variation => {
                // QUICKREF: All CIGAR diffs were 0 - output 0|0 immediately
                // Only use this path if --alignment-all is NOT set
                let mut repeat_mut = repeat.clone();
                if args.debug {
                    repeat_mut.set_time_stamp();
                }

                // Get reference sequence
                let fasta_reader = rust_htslib::faidx::Reader::from_path(&args.fasta)
                    .unwrap_or_else(|_| panic!("Failed to open FASTA file: {}", args.fasta));

                let repeat_ref_seq =
                    match repeat_mut.reference_repeat_sequence_with_reader(&fasta_reader) {
                        Some(seq) => seq,
                        None => {
                            error!("Cannot get reference sequence for {repeat}");
                            continue;
                        }
                    };

                let flags = vec!["QUICKREF".to_string()];
                let vcf_record = crate::vcf::VCFRecord::quick_reference(
                    &repeat_mut,
                    &repeat_ref_seq,
                    flags,
                    args,
                );
                results.push(vcf_record);
            }
            Some(info) => {
                // Needs full genotyping: extract sequences from filtered records
                let mut reads = parse_bam::Reads {
                    phase0: Vec::new(),
                    phase1: Vec::new(),
                    phase2: Vec::new(),
                    ps: None,
                };

                for &idx in &info.record_indices {
                    let record = &filtered_records[idx];
                    let phase = parse_bam::get_phase(record);
                    let seq = record.seq().as_bytes();

                    if args.unphased {
                        reads.phase0.push(seq);
                    } else {
                        match phase {
                            1 => {
                                reads.phase1.push(seq);
                                reads.ps = parse_bam::get_phase_set(record);
                            }
                            2 => {
                                reads.phase2.push(seq);
                                reads.ps = parse_bam::get_phase_set(record);
                            }
                            _ => {
                                reads.phase0.push(seq);
                            }
                        }
                    }
                }

                // Check if we have enough reads
                if reads.is_empty() {
                    if args.unphased {
                        debug!("Cannot genotype {repeat}: no reads found");
                    } else {
                        debug!(
                            "Cannot genotype {repeat}: no phased reads found. Use --unphased to genotype unphased reads."
                        );
                    }
                    // Use existing function to create missing genotype VCF record
                    let mut repeat_mut = repeat.clone();
                    if args.debug {
                        repeat_mut.set_time_stamp();
                    }
                    match genotype::genotype_repeat_singlethreaded(&mut repeat_mut, args, bam) {
                        Ok(vcf_record) => results.push(vcf_record),
                        Err(e) => {
                            error!("Problem processing {repeat}: {e}");
                        }
                    }
                    continue;
                }

                // Apply downsampling if needed
                if args.max_number_reads != -1 {
                    let max_reads = args.max_number_reads as usize;
                    parse_bam::downsample_reads_inplace(&mut reads.phase0, max_reads);
                    parse_bam::downsample_reads_inplace(&mut reads.phase1, max_reads / 2);
                    parse_bam::downsample_reads_inplace(&mut reads.phase2, max_reads / 2);
                }

                // Call the genotyping with the extracted reads
                let mut repeat_mut = repeat.clone();
                if args.debug {
                    repeat_mut.set_time_stamp();
                }
                match genotype::genotype_with_extracted_reads(&mut repeat_mut, reads, args) {
                    Ok(vcf_record) => results.push(vcf_record),
                    Err(e) => {
                        error!("Problem processing {repeat}: {e}");
                    }
                }
            }
        }
    }

    Ok(results)
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

    // Convert iterator to vector for batching
    let repeats_vec = repeats.collect();

    // Create batches of nearby targets (10KB default threshold)
    let batch_distance_threshold = 10000; // 10KB - targets within this distance are grouped
    let batches = create_batches(repeats_vec, batch_distance_threshold);

    info!("Processing {} targets in {} batches", num_intervals, batches.len());

    if args.threads == 1 {
        // When running single threaded things become easier and the tool will require less memory
        // Output is returned in the same order as the bed, and therefore not sorted before writing immediately to stdout
        // The indexedreader is created once and passed on to the function
        let mut bam = parse_bam::create_bam_reader(&args.bam, &args.fasta);
        for batch in &batches {
            match process_batch(batch, &mut bam, &args) {
                Ok(vcf_records) => {
                    for record in vcf_records {
                        writeln!(handle, "{record}").expect("Failed writing the result.");
                        pb.inc(1); // Increment progress for each locus
                    }
                }
                Err(e) => {
                    error!(
                        "Failed processing batch {}:{}-{}: {}",
                        batch.chromosome, batch.start, batch.end, e
                    );
                    // Still increment for failed loci
                    pb.inc(batch.repeats.len() as u64);
                }
            }
        }
    } else if args.sorted {
        // output is sorted by chrom, start and end
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .expect("Failed to create threadpool");

        // Create BAM reader pool for parallel processing
        let reader_pool = BamReaderPool::new(args.bam.clone(), args.fasta.clone());

        // Process all batches in parallel, collect results, then sort
        let genotypes = Mutex::new(Vec::with_capacity(num_intervals));

        batches.par_iter().for_each(|batch| {
            // Get thread-local BAM reader from pool
            reader_pool.with_reader(|reader| {
                match process_batch(batch, reader, &args) {
                    Ok(vcf_records) => {
                        let count = vcf_records.len() as u64;
                        let mut geno = genotypes.lock().expect("Unable to lock genotypes mutex");
                        geno.extend(vcf_records);
                        drop(geno); // Release lock before updating progress
                        pb.inc(count); // Increment progress for each locus
                    }
                    Err(e) => {
                        error!(
                            "Failed processing batch {}:{}-{}: {}",
                            batch.chromosome, batch.start, batch.end, e
                        );
                        pb.inc(batch.repeats.len() as u64);
                    }
                }
            });
        });

        let mut genotypes_vec = genotypes.lock().unwrap();
        // The final output is sorted by chrom, start and end
        genotypes_vec.sort_unstable();
        for g in &*genotypes_vec {
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

        // Create BAM reader pool for parallel processing
        let reader_pool = BamReaderPool::new(args.bam.clone(), args.fasta.clone());

        batches.par_iter().for_each(|batch| {
            // Get thread-local BAM reader from pool
            reader_pool.with_reader(|reader| {
                match process_batch(batch, reader, &args) {
                    Ok(vcf_records) => {
                        for record in vcf_records {
                            println!("{record}");
                            pb.inc(1); // Increment progress for each locus
                        }
                    }
                    Err(e) => {
                        error!(
                            "Failed processing batch {}:{}-{}: {}",
                            batch.chromosome, batch.start, batch.end, e
                        );
                        pb.inc(batch.repeats.len() as u64);
                    }
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
