use log::warn;
use rand::seq::SliceRandom;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::Aux;
use std::env;
use url::Url;

/// Net length difference from the reference over `[start - 1, end]`, counting every indel
/// and clip regardless of size. Used to decide whether a locus can be called homozygous
/// reference without aligning, where even a one-base difference disqualifies it.
///
/// `start` is 1-based (as in [`crate::repeats::RepeatInterval`]) while `reference_position`
/// walks in 0-based coordinates, so the bound is `start - 1`; the comparison is inclusive
/// at both ends because an aligner overwhelmingly places a repeat's indel on the repeat's
/// own first base, and excluding the boundary made the check blind to exactly the events it
/// exists to detect.
pub fn calculate_all_length_diff_from_cigar(
    record: &rust_htslib::bam::Record,
    start: u32,
    end: u32,
) -> i64 {
    let mut length_diff: i64 = 0;
    // the repeat occupies 0-based [start - 1, end); an indel on either boundary belongs to it
    let start = start.saturating_sub(1);
    let mut reference_position = record.reference_start() as u32;

    for entry in record.cigar().iter() {
        match entry {
            rust_htslib::bam::record::Cigar::Match(len)
            | rust_htslib::bam::record::Cigar::Equal(len)
            | rust_htslib::bam::record::Cigar::Diff(len) => {
                reference_position += *len;
            }
            rust_htslib::bam::record::Cigar::Del(len) => {
                if start <= reference_position && reference_position <= end {
                    length_diff -= i64::from(*len);
                }
                reference_position += *len;
            }
            rust_htslib::bam::record::Cigar::SoftClip(len)
                if start <= reference_position && reference_position <= end =>
            {
                length_diff += i64::from(*len);
            }
            rust_htslib::bam::record::Cigar::Ins(len)
                if start <= reference_position && reference_position <= end =>
            {
                length_diff += i64::from(*len);
            }
            rust_htslib::bam::record::Cigar::RefSkip(len) => reference_position += *len,
            _ => (),
        }
    }

    length_diff
}

/// Only insertions of at least this size may be taken in from the flank. One- and two-base
/// indels are everywhere in long-read alignments; counting them would let ordinary
/// sequencing noise beside the locus change the reported allele length. Weakly determined:
/// on the sample tested, anything from 3 to 10 performs within noise of the other.
const FLANK_MIN_INDEL: i64 = 3;
/// Minimum size of a clip that makes a read unusable near the locus.
const CLIP_MIN_LEN: i64 = 10;

/// A read whose alignment terminates (soft- or hard-clips) within this many bases of the
/// repeat cannot be trusted to have aligned through it: the clipped bases may be repeat
/// sequence the aligner gave up on, which is exactly what a large expansion looks like.
/// Such reads are dropped, which hands the locus to the alignment path once too few are
/// left - over-rejecting only costs time, under-rejecting costs a false reference call.
const CLIP_SEARCH: i64 = 100;

/// Slice the repeat allele straight out of the existing alignment.
///
/// Walks the CIGAR of `record` and returns the read bases the aligner placed over the
/// repeat interval: reference-matching repeat copies as well as insertions, shortened by
/// deletions. This is the same quantity the alignment-based path recovers as an insertion
/// against the repeat-compressed reference, but without re-aligning anything.
///
/// Aligners do not agree with each other, or between reads, about where a large insertion
/// belongs: at one locus the same ~900 bp insertion is placed anywhere across a 30 bp
/// stretch. Insertions of at least [`FLANK_MIN_INDEL`] bases lying within `flank` bases of
/// the interval are therefore counted towards the allele as well. The allele grows by
/// exactly their length; the bases kept are the ones nearest the repeat rather than the
/// inserted ones themselves, which leaves the length exact and the sequence approximate in
/// that (uncommon) case.
///
/// Deletions need no such tolerance: a deletion has a reference span, so one that belongs
/// to the repeat overlaps the interval and is already reflected in the slice. A deletion
/// lying entirely outside is a neighbouring event and must not shorten the allele.
///
/// Returns `None` when the read does not fully span the padded window, when it clips within
/// [`CLIP_SEARCH`] bases of the repeat, or when the repeat is deleted from the read
/// entirely; the alignment path likewise reports no insertion for such a read.
pub fn repeat_sequence_from_alignment(
    record: &bam::Record,
    start: u32,
    end: u32,
    flank: u32,
) -> Option<Vec<u8>> {
    use rust_htslib::bam::record::Cigar;

    // repeat.start is 1-based inclusive and repeat.end is the BED end, so the repeat
    // occupies 0-based [start - 1, end) - the same span the repeat-compressed reference
    // excises between its two flanks.
    let repeat_lo = i64::from(start) - 1;
    let repeat_hi = i64::from(end);
    let flank = i64::from(flank);
    // query offsets wanted for, in ascending order: the outer edge of the left flank, the
    // repeat start, the repeat end, and the outer edge of the right flank
    let targets = [repeat_lo - flank, repeat_lo, repeat_hi, repeat_hi + flank];
    let mut anchors = [None::<i64>; 4];
    // lengths of the qualifying insertions sitting in either flank
    let (mut left_insertions, mut right_insertions) = (0i64, 0i64);

    let mut ref_pos = record.reference_start();
    if ref_pos > targets[0] {
        return None; // read starts inside the window
    }
    let mut q_pos: i64 = 0;

    for entry in record.cigar().iter() {
        match entry {
            Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) => {
                let l = i64::from(*l);
                set_anchors(&mut anchors, &targets, ref_pos, l, q_pos, true);
                ref_pos += l;
                q_pos += l;
            }
            Cigar::Ins(l) => {
                let l = i64::from(*l);
                // an insertion on either boundary of the repeat already falls inside the
                // slice, so only the ones strictly out in a flank are counted here
                if l >= FLANK_MIN_INDEL {
                    if ref_pos >= targets[0] && ref_pos < targets[1] {
                        left_insertions += l;
                    } else if ref_pos > targets[2] && ref_pos <= targets[3] {
                        right_insertions += l;
                    }
                }
                q_pos += l;
            }
            Cigar::Del(l) | Cigar::RefSkip(l) => {
                let l = i64::from(*l);
                set_anchors(&mut anchors, &targets, ref_pos, l, q_pos, false);
                ref_pos += l;
            }
            Cigar::SoftClip(l) => {
                if i64::from(*l) >= CLIP_MIN_LEN && clips_into_locus(ref_pos, repeat_lo, repeat_hi)
                {
                    return None;
                }
                q_pos += i64::from(*l);
            }
            Cigar::HardClip(l) => {
                if i64::from(*l) >= CLIP_MIN_LEN && clips_into_locus(ref_pos, repeat_lo, repeat_hi)
                {
                    return None;
                }
            }
            Cigar::Pad(_) => (),
        }
        // stop once the window is behind us and no clip can still reach back to the locus
        if anchors[3].is_some() && ref_pos > repeat_hi + CLIP_SEARCH {
            break;
        }
    }
    // a read whose alignment ends flush with the window edge still spans it
    if anchors[3].is_none() && ref_pos >= targets[3] {
        anchors[3] = Some(q_pos);
    }

    let (win_lo, rep_lo, rep_hi, win_hi) = (anchors[0]?, anchors[1]?, anchors[2]?, anchors[3]?);
    // the allele is the repeat itself plus the insertions rescued from the flanks; take
    // those extra bases from the query sequence adjoining the repeat, without ever
    // reaching past the flank the aligner actually gave us
    let lo = (rep_lo - left_insertions).max(win_lo);
    let hi = (rep_hi + right_insertions).min(win_hi);
    if hi <= lo {
        return None;
    }
    let seq = record.seq();
    if hi as usize > seq.len() {
        return None;
    }
    // decode only the window: `as_bytes()` would unpack the whole (tens of kb) read for
    // the handful of bases the repeat actually needs
    Some(
        (lo as usize..hi as usize)
            .map(|i| seq[i].to_ascii_uppercase())
            .collect(),
    )
}

/// Record the query offset of every window edge this CIGAR operation reaches.
///
/// `q_pos` is the operation's query offset; an operation that consumes only reference (a
/// deletion) maps every reference position inside it to that same offset. The repeat start
/// and the left window edge are claimed as soon as the operation reaches them, so that an
/// insertion sitting on that boundary falls inside the slice; the repeat end and the right
/// edge are claimed only once the operation moves past them, for the same reason.
fn set_anchors(
    anchors: &mut [Option<i64>; 4],
    targets: &[i64; 4],
    ref_pos: i64,
    len: i64,
    q_pos: i64,
    consumes_query: bool,
) {
    for (i, (anchor, target)) in anchors.iter_mut().zip(targets.iter()).enumerate() {
        let reached = if i < 2 {
            ref_pos + len >= *target
        } else {
            ref_pos + len > *target
        };
        if anchor.is_none() && reached {
            *anchor = Some(if consumes_query {
                q_pos + (*target - ref_pos).max(0)
            } else {
                q_pos
            });
        }
    }
}

/// Whether an alignment that terminates at `ref_pos` does so close enough to the repeat
/// that its clipped bases might be repeat sequence.
fn clips_into_locus(ref_pos: i64, repeat_lo: i64, repeat_hi: i64) -> bool {
    ref_pos >= repeat_lo - CLIP_SEARCH && ref_pos <= repeat_hi + CLIP_SEARCH
}

pub struct Reads {
    pub phase0: Vec<Vec<u8>>, // Unphased or haploid reads
    pub phase1: Vec<Vec<u8>>, // Phase 1 reads (phased mode)
    pub phase2: Vec<Vec<u8>>, // Phase 2 reads (phased mode)
    pub ps: Option<u32>,
    /// True when the stored sequences are already the repeat alleles sliced out of the
    /// alignment (`--fast`), rather than whole reads still to be re-aligned.
    pub sliced: bool,
}

impl Reads {
    /// Returns true if all phase vectors are empty
    pub fn is_empty(&self) -> bool {
        self.phase0.is_empty() && self.phase1.is_empty() && self.phase2.is_empty()
    }
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

/// Downsample a vector of reads to a maximum number in-place
pub fn downsample_reads_inplace(phase_reads: &mut Vec<Vec<u8>>, max_reads: usize) {
    let n_reads = phase_reads.len();
    if n_reads <= max_reads {
        return; // No downsampling needed
    }

    // Use partial_shuffle for efficient random selection
    // Shuffles the first max_reads elements randomly from the full vector
    let mut rng = rand::rng();
    phase_reads.partial_shuffle(&mut rng, max_reads);

    // Keep only the first max_reads elements
    phase_reads.truncate(max_reads);
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
    // Create struct and populate directly
    let mut reads = Reads {
        phase0: Vec::new(),
        phase1: Vec::new(),
        phase2: Vec::new(),
        ps: None,
        sliced: false,
    };

    // Extract sequences from collected records
    for r in records {
        if unphased {
            // if unphased put reads in phase 0
            reads.phase0.push(r.seq().as_bytes());
        } else {
            let phase = get_phase(&r);
            match phase {
                1 => {
                    reads.phase1.push(r.seq().as_bytes());
                }
                2 => {
                    reads.phase2.push(r.seq().as_bytes());
                }
                _ => {} // phase 0 or invalid - skip
            }
            reads.ps = get_phase_set(&r);
        }
    }

    if reads.is_empty() {
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
        Some(reads)
    } else {
        // if more than <max_number_reads> spanning reads are found, randomly select just <max_number_reads> items
        // if unphased, just select <max_number_reads> reads from phase 0
        // if phased, select <max_number_reads>/2 reads from each phase
        let max_number_reads = max_number_reads as usize;

        downsample_reads_inplace(&mut reads.phase0, max_number_reads);
        downsample_reads_inplace(&mut reads.phase1, max_number_reads / 2);
        downsample_reads_inplace(&mut reads.phase2, max_number_reads / 2);

        Some(reads)
    }
}

pub fn get_overlapping_reads(
    bam: &mut bam::IndexedReader,
    repeat: &crate::repeats::RepeatInterval,
    unphased: bool,
    max_number_reads: isize,
    min_mapq: u8,
) -> Option<Reads> {
    let tid = bam
        .header()
        .tid(repeat.chrom.as_bytes())
        .unwrap_or_else(|| panic!("Invalid chromosome {}", repeat.chrom));
    // repeat.start is 1-based, but bam.fetch expects 0-based half-open coordinates
    bam.fetch((tid, repeat.start - 1, repeat.end))
        .unwrap_or_else(|err| panic!("Failure to extract reads from bam for {repeat}:\n{err}"));

    // Create struct and populate directly
    let mut reads = Reads {
        phase0: Vec::new(),
        phase1: Vec::new(),
        phase2: Vec::new(),
        ps: None,
        sliced: false,
    };

    // extract sequences spanning the repeat locus
    for r in bam.rc_records() {
        let r = r.unwrap_or_else(|err| panic!("Error reading BAM file in region {repeat}:\n{err}"));
        // skip poorly mapped reads and reads that do not span the repeat locus
        if r.mapq() < min_mapq
            || r.reference_start() > repeat.start.into()
            || r.reference_end() < repeat.end.into()
        {
            continue;
        }
        if unphased {
            // if unphased put reads in phase 0
            reads.phase0.push(r.seq().as_bytes());
        } else {
            let phase = get_phase(&r);
            match phase {
                1 => {
                    reads.phase1.push(r.seq().as_bytes());
                    reads.ps = get_phase_set(&r);
                }
                2 => {
                    reads.phase2.push(r.seq().as_bytes());
                    reads.ps = get_phase_set(&r);
                }
                _ => {} // phase 0 or invalid - skip
            }
        }
    }

    if reads.is_empty() {
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
        Some(reads)
    } else {
        // if more than <max_number_reads> spanning reads are found, randomly select just <max_number_reads> items
        // if unphased, just select <max_number_reads> reads from phase 0
        // if phased, select <max_number_reads>/2 reads from each phase
        let max_number_reads = max_number_reads as usize;

        downsample_reads_inplace(&mut reads.phase0, max_number_reads);
        downsample_reads_inplace(&mut reads.phase1, max_number_reads / 2);
        downsample_reads_inplace(&mut reads.phase2, max_number_reads / 2);

        Some(reads)
    }
}

pub fn get_phase(record: &bam::Record) -> u8 {
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

pub fn get_phase_set(record: &bam::Record) -> Option<u32> {
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
    let _reads = get_overlapping_reads(&mut bam, &repeat, unphased, 60, 10);
}

#[test]
#[ignore = "requires network access to ftp.1000genomes.ebi.ac.uk - set TEST_REMOTE_BAM=1 to enable"]
fn test_get_overlapping_reads_url() {
    if std::env::var("TEST_REMOTE_BAM").is_err() {
        return;
    }
    let bam = String::from(
        "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/hg38/HG00096.hg38.cram",
    );
    let fasta = String::from("test_data/chr7.fa.gz");
    let repeat = crate::repeats::RepeatInterval {
        chrom: String::from("chr7"),
        start: 154654404,
        end: 154654432,
        created: None,
    };
    // this sample is aligned without HP tags, so only the unphased path returns reads
    let unphased = true;
    let mut bam = create_bam_reader(&bam, &fasta);
    let reads = get_overlapping_reads(&mut bam, &repeat, unphased, 60, 10)
        .expect("Expected spanning reads from the remote CRAM");
    // the locus had 20 spanning MAPQ60 reads when this test was written; assert loosely
    // so a re-release of the file upstream does not break the test
    assert!(
        reads.phase0.len() >= 10,
        "Expected at least 10 spanning reads from the remote CRAM, got {}",
        reads.phase0.len()
    );
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

#[cfg(test)]
mod length_diff_tests {
    use super::*;
    use rust_htslib::bam::record::{Cigar, CigarString, Record};

    fn record(pos: i64, cigar: Vec<Cigar>) -> Record {
        let mut rec = Record::new();
        let query_len: u32 = cigar
            .iter()
            .map(|op| match op {
                Cigar::Match(l)
                | Cigar::Equal(l)
                | Cigar::Diff(l)
                | Cigar::Ins(l)
                | Cigar::SoftClip(l) => *l,
                _ => 0,
            })
            .sum();
        let seq = vec![b'A'; query_len as usize];
        let qual = vec![30u8; query_len as usize];
        rec.set(b"read", Some(&CigarString(cigar)), &seq, &qual);
        rec.set_pos(pos);
        rec
    }

    // repeat at 1-based 1001-1010, i.e. 0-based [1000, 1010)
    const START: u32 = 1001;
    const END: u32 = 1010;

    #[test]
    fn counts_an_insertion_at_the_first_base_of_the_repeat() {
        // Aligners overwhelmingly park a repeat's insertion on the repeat's first base.
        // Missing it here made QUICKREF call such loci homozygous reference.
        let rec = record(950, vec![Cigar::Match(50), Cigar::Ins(6), Cigar::Match(50)]);
        assert_eq!(calculate_all_length_diff_from_cigar(&rec, START, END), 6);
    }

    #[test]
    fn counts_an_indel_at_the_last_base_of_the_repeat() {
        let rec = record(950, vec![Cigar::Match(59), Cigar::Ins(4), Cigar::Match(41)]);
        assert_eq!(calculate_all_length_diff_from_cigar(&rec, START, END), 4);
        let rec = record(950, vec![Cigar::Match(59), Cigar::Del(4), Cigar::Match(41)]);
        assert_eq!(calculate_all_length_diff_from_cigar(&rec, START, END), -4);
    }

    #[test]
    fn counts_an_indel_in_the_middle_of_the_repeat() {
        let rec = record(950, vec![Cigar::Match(55), Cigar::Del(3), Cigar::Match(45)]);
        assert_eq!(calculate_all_length_diff_from_cigar(&rec, START, END), -3);
    }

    #[test]
    fn ignores_an_indel_outside_the_window() {
        let rec = record(950, vec![Cigar::Match(20), Cigar::Ins(9), Cigar::Match(80)]);
        assert_eq!(calculate_all_length_diff_from_cigar(&rec, START, END), 0);
    }

    #[test]
    fn a_clean_read_reports_no_difference() {
        let rec = record(950, vec![Cigar::Match(100)]);
        assert_eq!(calculate_all_length_diff_from_cigar(&rec, START, END), 0);
    }
}

#[cfg(test)]
mod slice_tests {
    use super::*;
    use rust_htslib::bam::record::{Cigar, CigarString, Record};

    /// Build a record at 0-based `pos` with `cigar` over `seq`.
    fn record(pos: i64, cigar: Vec<Cigar>, seq: &[u8]) -> Record {
        let mut rec = Record::new();
        let qual = vec![30u8; seq.len()];
        rec.set(b"read", Some(&CigarString(cigar)), seq, &qual);
        rec.set_pos(pos);
        rec
    }

    /// A repeat at 1-based start=1001, end=1010 occupies 0-based [1000, 1010).
    const START: u32 = 1001;
    const END: u32 = 1010;

    #[test]
    fn slices_reference_matching_repeat() {
        // 100M from 0-based 950: the repeat window is plain reference sequence
        let mut seq = vec![b'A'; 100];
        seq[50..60].copy_from_slice(b"CGCGCGCGCG"); // 0-based 1000..1010
        let rec = record(950, vec![Cigar::Match(100)], &seq);
        let out = repeat_sequence_from_alignment(&rec, START, END, 0).unwrap();
        assert_eq!(out, b"CGCGCGCGCG".to_vec());
    }

    #[test]
    fn includes_insertion_inside_the_repeat() {
        // 55M 6I 45M from 950: a 6bp insertion at 0-based 1005, inside the repeat
        let mut seq = vec![b'A'; 106];
        seq[50..55].copy_from_slice(b"CGCGC");
        seq[55..61].copy_from_slice(b"TTTTTT");
        seq[61..66].copy_from_slice(b"GCGCG");
        let rec = record(950, vec![Cigar::Match(55), Cigar::Ins(6), Cigar::Match(45)], &seq);
        let out = repeat_sequence_from_alignment(&rec, START, END, 0).unwrap();
        assert_eq!(out, b"CGCGCTTTTTTGCGCG".to_vec());
    }

    #[test]
    fn deletion_shortens_the_allele() {
        // 55M 4D 45M from 950: 4 reference bases of the repeat are missing from the read
        let mut seq = vec![b'A'; 100];
        seq[50..55].copy_from_slice(b"CGCGC");
        seq[55..61].copy_from_slice(b"GCGCGC");
        let rec = record(950, vec![Cigar::Match(55), Cigar::Del(4), Cigar::Match(45)], &seq);
        let out = repeat_sequence_from_alignment(&rec, START, END, 0).unwrap();
        assert_eq!(out, b"CGCGCG".to_vec(), "10bp repeat minus a 4bp deletion");
    }

    #[test]
    fn flank_captures_a_large_insertion_shifted_out_of_the_interval() {
        // 12bp insertion at 0-based 996, four bases before the annotated repeat start
        let mut seq = vec![b'A'; 112];
        seq[46..58].copy_from_slice(b"TTTTTTTTTTTT"); // the insertion
        seq[62..72].copy_from_slice(b"CGCGCGCGCG"); // 0-based 1000..1010
        let rec = record(950, vec![Cigar::Match(46), Cigar::Ins(12), Cigar::Match(54)], &seq);
        // without flank only the annotated interval is taken
        let tight = repeat_sequence_from_alignment(&rec, START, END, 0).unwrap();
        assert_eq!(tight, b"CGCGCGCGCG".to_vec());
        // with flank the allele grows by exactly the 12 inserted bases
        let padded = repeat_sequence_from_alignment(&rec, START, END, 30).unwrap();
        assert_eq!(padded.len(), 22);
        assert!(padded.ends_with(b"CGCGCGCGCG"));
    }

    #[test]
    fn small_flanking_insertion_is_noise_and_is_ignored() {
        // 1bp insertion at 0-based 1015, five bases past the repeat: ordinary ONT error
        let mut seq = vec![b'A'; 101];
        seq[50..60].copy_from_slice(b"CGCGCGCGCG");
        let rec = record(950, vec![Cigar::Match(65), Cigar::Ins(1), Cigar::Match(35)], &seq);
        for flank in [0, 10, 30] {
            let out = repeat_sequence_from_alignment(&rec, START, END, flank).unwrap();
            assert_eq!(out, b"CGCGCGCGCG".to_vec(), "flank {flank}");
        }
    }

    #[test]
    fn deletion_outside_the_interval_does_not_shorten_the_allele() {
        // 5bp deletion at 0-based 1012..1017, entirely past the repeat. It belongs to the
        // flank, not to the repeat, so the allele must stay at reference length.
        let mut seq = vec![b'A'; 95];
        seq[50..60].copy_from_slice(b"CGCGCGCGCG");
        let rec = record(950, vec![Cigar::Match(62), Cigar::Del(5), Cigar::Match(33)], &seq);
        for flank in [0, 10, 30] {
            let out = repeat_sequence_from_alignment(&rec, START, END, flank).unwrap();
            assert_eq!(out, b"CGCGCGCGCG".to_vec(), "flank {flank}");
        }
    }

    #[test]
    fn deletion_outside_the_interval_never_discards_the_read() {
        // a flanking deletion longer than the repeat itself must not empty the slice
        let mut seq = vec![b'A'; 88];
        seq[50..60].copy_from_slice(b"CGCGCGCGCG");
        let rec = record(950, vec![Cigar::Match(62), Cigar::Del(12), Cigar::Match(26)], &seq);
        let out = repeat_sequence_from_alignment(&rec, START, END, 30);
        assert_eq!(out, Some(b"CGCGCGCGCG".to_vec()));
    }

    #[test]
    fn rejects_a_read_clipped_near_but_outside_the_window() {
        // 900bp soft clip ending 11 bases before the repeat: the clipped bases may well be
        // the expansion the aligner refused to align through, so the read is unusable
        let seq = vec![b'A'; 1000];
        let rec = record(989, vec![Cigar::SoftClip(900), Cigar::Match(100)], &seq);
        assert!(repeat_sequence_from_alignment(&rec, START, END, 10).is_none());
        // the same read hard-clipped, as a supplementary alignment would carry it
        let rec = record(989, vec![Cigar::HardClip(900), Cigar::Match(100)], &seq[..100]);
        assert!(repeat_sequence_from_alignment(&rec, START, END, 10).is_none());
        // a short clip far from the locus is harmless
        let mut seq = vec![b'A'; 1000];
        seq[50..60].copy_from_slice(b"CGCGCGCGCG");
        let rec = record(950, vec![Cigar::Match(500), Cigar::SoftClip(500)], &seq);
        assert_eq!(
            repeat_sequence_from_alignment(&rec, START, END, 10),
            Some(b"CGCGCGCGCG".to_vec())
        );
    }

    #[test]
    fn rejects_reads_that_do_not_span() {
        let seq = vec![b'A'; 100];
        // read starts inside the padded window
        let rec = record(1005, vec![Cigar::Match(100)], &seq);
        assert!(repeat_sequence_from_alignment(&rec, START, END, 0).is_none());
        // read ends inside the window
        let rec = record(950, vec![Cigar::Match(55)], &seq[..55]);
        assert!(repeat_sequence_from_alignment(&rec, START, END, 0).is_none());
        // read is soft-clipped inside the window
        let rec = record(950, vec![Cigar::Match(55), Cigar::SoftClip(45)], &seq);
        assert!(repeat_sequence_from_alignment(&rec, START, END, 0).is_none());
    }
}
