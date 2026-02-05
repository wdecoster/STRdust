use log::warn;
use rand::seq::IndexedRandom;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::Aux;
use std::collections::HashMap;
use std::env;
use url::Url;

/// Calculate length differences from reference for quick reference check
/// Counts ALL indels/clips regardless of size within the specified region
/// Used for fast 0|0 detection where even small variants matter
pub fn calculate_all_length_diff_from_cigar(
    record: &rust_htslib::bam::Record,
    start: u32,
    end: u32,
) -> i64 {
    let mut length_diff: i64 = 0;
    let mut reference_position = record.reference_start() as u32;

    for entry in record.cigar().iter() {
        match entry {
            rust_htslib::bam::record::Cigar::Match(len)
            | rust_htslib::bam::record::Cigar::Equal(len)
            | rust_htslib::bam::record::Cigar::Diff(len) => {
                reference_position += *len;
            }
            rust_htslib::bam::record::Cigar::Del(len) => {
                if start < reference_position && reference_position < end {
                    length_diff -= i64::from(*len);
                }
                reference_position += *len;
            }
            rust_htslib::bam::record::Cigar::SoftClip(len) => {
                if start < reference_position && reference_position < end {
                    length_diff += i64::from(*len);
                }
            }
            rust_htslib::bam::record::Cigar::Ins(len) => {
                if start < reference_position && reference_position < end {
                    length_diff += i64::from(*len);
                }
            }
            rust_htslib::bam::record::Cigar::RefSkip(len) => reference_position += *len,
            _ => (),
        }
    }

    length_diff
}

pub struct Reads {
    // could consider not to use a hashmap here and use an attribute per phase
    pub seqs: HashMap<u8, Vec<Vec<u8>>>,
    pub ps: Option<u32>,
}

/// Sets up the CURL_CA_BUNDLE environment variable for HTTPS/S3 access
/// Tries to use a CA bundle from standard locations, with appropriate fallbacks
fn setup_ssl_certificates() {
    // Only configure if not already set by the user
    if env::var("CURL_CA_BUNDLE").is_ok() {
        return;
    }

    // Common CA bundle locations across different systems
    let possible_paths = vec![
        "/etc/ssl/certs/ca-certificates.crt",     // Debian/Ubuntu
        "/etc/pki/tls/certs/ca-bundle.crt",       // RHEL/CentOS/Amazon Linux
        "/etc/ssl/ca-bundle.pem",                 // SUSE
        "/usr/local/share/certs/ca-root-nss.crt", // FreeBSD
        "/usr/local/etc/openssl/cert.pem",        // macOS Homebrew
        "/etc/ssl/cert.pem",                      // macOS/OpenBSD
    ];

    // Try each path in order
    for path in possible_paths {
        if std::path::Path::new(path).exists() {
            unsafe { env::set_var("CURL_CA_BUNDLE", path) };
            return;
        }
    }

    // None of the paths exist, warn the user
    warn!(
        "Could not find a valid CA certificate bundle for HTTPS/S3 access. \
        HTTPS/S3 connections may fail. Set the CURL_CA_BUNDLE environment \
        variable to the path of your system's CA certificate bundle."
    );
}

pub fn create_bam_reader(bamf: &str, fasta: &str) -> bam::IndexedReader {
    let mut bam = if bamf.starts_with("s3") || bamf.starts_with("https://") {
        setup_ssl_certificates();
        bam::IndexedReader::from_url(&Url::parse(bamf).expect("Failed to parse URL"))
            .unwrap_or_else(|err| panic!("Error opening remote BAM: {err}"))
    } else {
        bam::IndexedReader::from_path(bamf)
            .unwrap_or_else(|err| panic!("Error opening local BAM: {err}"))
    };
    if bamf.ends_with(".cram") {
        bam.set_cram_options(
            hts_sys::hts_fmt_option_CRAM_OPT_REQUIRED_FIELDS,
            hts_sys::sam_fields_SAM_AUX
                | hts_sys::sam_fields_SAM_MAPQ
                | hts_sys::sam_fields_SAM_SEQ,
        )
        .expect("Failed setting cram options");
        bam.set_reference(fasta)
            .expect("Failed setting reference for CRAM file");
    }
    bam
}

/// Process already-collected BAM records into phased read sequences
///
/// This function takes a vector of pre-fetched BAM records and extracts their sequences
/// according to phase, avoiding a second BAM fetch operation. Used by the fast reference
/// check path to eliminate redundant I/O.
///
/// # Arguments
/// * `records` - Pre-collected BAM records
/// * `repeat` - STR target region (for logging)
/// * `unphased` - Whether to treat reads as unphased
/// * `max_number_reads` - Maximum reads per haplotype (-1 = all reads)
///
/// # Returns
/// Some(Reads) with phased sequences, or None if no reads found
pub fn process_collected_reads(
    records: Vec<bam::Record>,
    repeat: &crate::repeats::RepeatInterval,
    unphased: bool,
    max_number_reads: isize,
) -> Option<Reads> {
    // Per haplotype the read sequences are kept in a dictionary
    let mut seqs = HashMap::from([(0, Vec::new()), (1, Vec::new()), (2, Vec::new())]);
    let mut ps = None;

    // Extract sequences from collected records
    for r in records {
        if unphased {
            // if unphased put reads in phase 0
            seqs.get_mut(&0).unwrap().push(r.seq().as_bytes());
        } else {
            let phase = get_phase(&r);
            if phase > 0 {
                let seq = r.seq().as_bytes();
                seqs.get_mut(&phase).unwrap().push(seq);
                ps = get_phase_set(&r);
            }
        }
    }

    if seqs.is_empty() {
        // error/warning message depends on whether we are looking for phased reads or not
        if unphased {
            warn!("Cannot genotype {repeat}: no reads found");
        } else {
            warn!(
                "Cannot genotype {repeat}: no phased reads found. Use --unphased to genotype unphased reads."
            );
        }
        None
    } else if max_number_reads == -1 {
        // if max_number_reads is -1, use all reads without downsampling
        Some(Reads { seqs, ps })
    } else {
        // if more than <max_number_reads> spanning reads are found in seqs, randomly select just <max_number_reads> items from the vector
        // if unphased, just select <max_number_reads> reads from phase 0
        // if phased, select <max_number_reads>/2 reads from each phase
        let max_number_reads = max_number_reads as usize;
        let mut rng = rand::rng();

        // Pre-allocate filtered vectors with exact size needed
        let mut seqs_filtered = HashMap::from([
            (0, Vec::with_capacity(max_number_reads)),
            (1, Vec::with_capacity(max_number_reads / 2)),
            (2, Vec::with_capacity(max_number_reads / 2)),
        ]);
        let max_reads_per_phase = HashMap::from([
            (0, max_number_reads),
            (1, max_number_reads / 2),
            (2, max_number_reads / 2),
        ]);
        for (phase, seqs_phase) in seqs.iter() {
            let n_reads = seqs_phase.len();
            let n_reads_to_select = n_reads.min(max_reads_per_phase[phase]);

            if n_reads_to_select > 0 {
                let selected_reads = seqs_phase.choose_multiple(&mut rng, n_reads_to_select);
                // Reserve exact capacity to avoid reallocation
                let phase_vec = seqs_filtered.get_mut(phase).unwrap();
                phase_vec.reserve_exact(n_reads_to_select);
                for read in selected_reads {
                    phase_vec.push(read.to_vec());
                }
            }
        }

        Some(Reads { seqs: seqs_filtered, ps })
    }
}

pub fn get_overlapping_reads(
    bam: &mut bam::IndexedReader,
    repeat: &crate::repeats::RepeatInterval,
    unphased: bool,
    max_number_reads: isize,
) -> Option<Reads> {
    let tid = bam
        .header()
        .tid(repeat.chrom.as_bytes())
        .unwrap_or_else(|| panic!("Invalid chromosome {}", repeat.chrom));
    // repeat.start is 1-based, but bam.fetch expects 0-based half-open coordinates
    bam.fetch((tid, repeat.start - 1, repeat.end))
        .unwrap_or_else(|err| panic!("Failure to extract reads from bam for {repeat}:\n{err}"));
    // Per haplotype the read sequences are kept in a dictionary
    let mut seqs = HashMap::from([(0, Vec::new()), (1, Vec::new()), (2, Vec::new())]);
    let mut ps = None;
    // extract sequences spanning the repeat locus
    for r in bam.rc_records() {
        let r = r.unwrap_or_else(|err| panic!("Error reading BAM file in region {repeat}:\n{err}"));
        // skip reads with mapq 0 or reads that do not span the repeat locus
        if r.mapq() == 0
            || r.reference_start() > repeat.start.into()
            || r.reference_end() < repeat.end.into()
        {
            continue;
        }
        if unphased {
            // if unphased put reads in phase 0
            seqs.get_mut(&0).unwrap().push(r.seq().as_bytes());
        } else {
            let phase = get_phase(&r);
            if phase > 0 {
                let seq = r.seq().as_bytes();
                seqs.get_mut(&phase).unwrap().push(seq);
                ps = get_phase_set(&r);
                // writing fasta to stdout
                // println!(">read_{}\n{}", phase, std::str::from_utf8(&seq).unwrap());
            }
        }
    }
    if seqs.is_empty() {
        // error/warning message depends on whether we are looking for phased reads or not
        if unphased {
            warn!("Cannot genotype {repeat}: no reads found");
        } else {
            warn!(
                "Cannot genotype {repeat}: no phased reads found. Use --unphased to genotype unphased reads."
            );
        }
        None
    } else if max_number_reads == -1 {
        // if max_number_reads is -1, use all reads without downsampling
        Some(Reads { seqs, ps })
    } else {
        // if more than <max_number_reads> spanning reads are found in seqs, randomly select just <max_number_reads> items from the vector
        // if unphased, just select <max_number_reads> reads from phase 0
        // if phased, select <max_number_reads>/2 reads from each phase
        let max_number_reads = max_number_reads as usize;
        let mut rng = rand::rng();

        // Pre-allocate filtered vectors with exact size needed
        let mut seqs_filtered = HashMap::from([
            (0, Vec::with_capacity(max_number_reads)),
            (1, Vec::with_capacity(max_number_reads / 2)),
            (2, Vec::with_capacity(max_number_reads / 2)),
        ]);
        let max_reads_per_phase = HashMap::from([
            (0, max_number_reads),
            (1, max_number_reads / 2),
            (2, max_number_reads / 2),
        ]);
        for (phase, seqs_phase) in seqs.iter() {
            let n_reads = seqs_phase.len();
            let n_reads_to_select = n_reads.min(max_reads_per_phase[phase]);

            if n_reads_to_select > 0 {
                let selected_reads = seqs_phase.choose_multiple(&mut rng, n_reads_to_select);
                // Reserve exact capacity to avoid reallocation
                let phase_vec = seqs_filtered.get_mut(phase).unwrap();
                phase_vec.reserve_exact(n_reads_to_select);
                for read in selected_reads {
                    phase_vec.push(read.to_vec());
                }
            }
        }

        Some(Reads { seqs: seqs_filtered, ps })
    }
}

fn get_phase(record: &bam::Record) -> u8 {
    match record.aux(b"HP") {
        Ok(value) => match value {
            Aux::U8(v) => v,
            Aux::U16(v) => u8::try_from(v).expect("Unexpected phase identifier for HP: {v:?}"),
            Aux::I32(v) => u8::try_from(v).expect("Unexpected phase identifier for HP: {v:?}"),
            _ => panic!("Unexpected type of Aux {value:?} for HP"),
        },
        Err(_e) => 0,
    }
}

fn get_phase_set(record: &bam::Record) -> Option<u32> {
    match record.aux(b"PS") {
        Ok(value) => match value {
            Aux::U32(v) => Some(v),
            Aux::I8(v) => {
                Some(u32::try_from(v).expect("Unexpected phase set identifier for PS: {v:?}"))
            }
            Aux::I16(v) => {
                Some(u32::try_from(v).expect("Unexpected phase set identifier for PS: {v:?}"))
            }
            Aux::I32(v) => {
                Some(u32::try_from(v).expect("Unexpected phase set identifier for PS: {v:?}"))
            }
            Aux::U8(v) => Some(u32::from(v)),
            Aux::U16(v) => Some(u32::from(v)),
            _ => panic!("Unexpected type of Aux {value:?} for PS"),
        },
        Err(_e) => None,
    }
}

#[test]
fn test_get_overlapping_reads() {
    let bam = String::from("test_data/small-test-phased.bam");
    let fasta = String::from("test_data/chr7.fa.gz");
    let repeat = crate::repeats::RepeatInterval {
        chrom: String::from("chr7"),
        start: 154654404,
        end: 154654432,
        created: None,
    };
    let unphased = false;
    let mut bam = create_bam_reader(&bam, &fasta);
    let _reads = get_overlapping_reads(&mut bam, &repeat, unphased, 60);
}

#[test]
fn test_get_overlapping_reads_url() {
    let bam = String::from(
        "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/hg38/HG00096.hg38.cram",
    );
    let fasta = String::from("test_data/chr7.fa.gz");
    let repeat = crate::repeats::RepeatInterval {
        chrom: String::from("chr20"),
        start: 154654404,
        end: 154654432,
        created: None,
    };
    let unphased = false;
    let mut bam = create_bam_reader(&bam, &fasta);
    let _reads = get_overlapping_reads(&mut bam, &repeat, unphased, 60);
}

#[test]
fn test_get_phase() {
    let mut bam =
        bam::Reader::from_path("test_data/small-test-phased.bam").expect("Failed opening bam");
    let record = bam
        .records()
        .next()
        .expect("Failed to read first record from bam");
    let phase = get_phase(&record.unwrap());
    assert_eq!(phase, 2);
}
