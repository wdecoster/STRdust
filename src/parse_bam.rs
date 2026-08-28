use log::warn;
use rand::seq::SliceRandom;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::Aux;
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
            rust_htslib::bam::record::Cigar::SoftClip(len)
                if start < reference_position && reference_position < end =>
            {
                length_diff += i64::from(*len);
            }
            rust_htslib::bam::record::Cigar::Ins(len)
                if start < reference_position && reference_position < end =>
            {
                length_diff += i64::from(*len);
            }
            rust_htslib::bam::record::Cigar::RefSkip(len) => reference_position += *len,
            _ => (),
        }
    }

    length_diff
}

/// How far beyond the annotated interval an indel is still considered part of the repeat,
/// mirroring the +/-30 base window the alignment path allows around its junction.
const FLANK_SEARCH: i64 = 30;
/// Only indels of at least this size widen the flank. Single-base indels are everywhere in
/// long-read alignments and would drag unrelated sequencing noise into the allele.
const FLANK_MIN_INDEL: i64 = 10;

/// Pick the flank to use on either side of the repeat for one read.
///
/// Starts from the requested `flank` and widens it just far enough to take in whole indels
/// of at least [`FLANK_MIN_INDEL`] bases lying within [`FLANK_SEARCH`] of the interval, so
/// that a large indel the aligner placed a little off the annotated boundary still counts
/// towards the allele.
fn flanks_for_read(record: &bam::Record, repeat_lo: i64, repeat_hi: i64, flank: i64) -> (i64, i64) {
    use rust_htslib::bam::record::Cigar;

    let (mut flank_lo, mut flank_hi) = (flank, flank);
    let search_lo = repeat_lo - FLANK_SEARCH;
    let search_hi = repeat_hi + FLANK_SEARCH;
    let mut ref_pos = record.reference_start();

    for entry in record.cigar().iter() {
        let (len, consumes_ref) = match entry {
            Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) => (i64::from(*l), true),
            Cigar::Del(l) | Cigar::RefSkip(l) => (i64::from(*l), true),
            Cigar::Ins(l) => (i64::from(*l), false),
            Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => (0, false),
        };
        let is_indel = matches!(entry, Cigar::Ins(_) | Cigar::Del(_));
        if is_indel && len >= FLANK_MIN_INDEL {
            // reference span of the event: insertions occupy no reference, deletions do
            let (event_lo, event_hi) = if consumes_ref {
                (ref_pos, ref_pos + len)
            } else {
                (ref_pos, ref_pos)
            };
            if event_lo >= search_lo && event_lo < repeat_lo {
                flank_lo = flank_lo.max(repeat_lo - event_lo);
            }
            if event_hi > repeat_hi && event_hi <= search_hi {
                flank_hi = flank_hi.max(event_hi - repeat_hi);
            }
        }
        if consumes_ref {
            ref_pos += len;
            if ref_pos > search_hi {
                break;
            }
        }
    }
    (flank_lo, flank_hi)
}

/// Slice the repeat allele straight out of the existing alignment.
///
/// Walks the CIGAR of `record` and returns the read bases the aligner placed over the
/// repeat interval: reference-matching repeat copies as well as insertions, shortened by
/// deletions. This is the same quantity the alignment-based path recovers as an insertion
/// against the repeat-compressed reference, but without re-aligning anything.
///
/// `flank` bases of reference on either side are taken into the window and then trimmed
/// back off again, so the allele grows or shrinks by exactly the indels the aligner placed
/// just outside the annotated interval - repeat boundaries in a catalog rarely agree with
/// where an aligner puts an indel. The trim removes `flank` *reference* bases, so a
/// flanking insertion still lengthens the allele; the bases kept are then the ones nearest
/// the repeat rather than the inserted ones themselves, but the allele length is exact
/// either way. With `flank = 0` the slice is exactly the annotated interval.
///
/// Aligners do not agree with each other, or between reads, about where a large indel
/// belongs: at one locus below the same ~900 bp insertion is placed anywhere across a 30 bp
/// stretch. The flank is therefore widened per read, up to [`FLANK_SEARCH`], to take in
/// whole indels of at least [`FLANK_MIN_INDEL`] bases that sit just outside the interval -
/// the same tolerance the alignment path applies around its junction.
///
/// Returns `None` when the read does not fully span the padded window (it starts or ends
/// inside it, or is soft/hard-clipped into it): such a read cannot report an allele length.
/// Also `None` when nothing is left after trimming, i.e. the read has deleted the whole
/// repeat; the alignment-based path likewise reports no insertion for such a read.
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
    let (flank_lo, flank_hi) = flanks_for_read(record, repeat_lo, repeat_hi, i64::from(flank));
    // query offsets wanted for the two ends of the padded window
    let targets = [repeat_lo - flank_lo, repeat_hi + flank_hi];
    let mut anchors = [None::<i64>; 2];

    let mut ref_pos = record.reference_start();
    if ref_pos > targets[0] {
        return None; // read starts inside the window
    }
    let mut q_pos: i64 = 0;

    for entry in record.cigar().iter() {
        match entry {
            Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) => {
                let l = i64::from(*l);
                // the lower anchor is placed as soon as the window is reached, so an
                // insertion sitting on that boundary falls inside the slice; the upper
                // anchor waits until the window is left behind, for the same reason
                if anchors[0].is_none() && ref_pos + l >= targets[0] {
                    anchors[0] = Some(q_pos + (targets[0] - ref_pos).max(0));
                }
                if anchors[1].is_none() && ref_pos + l > targets[1] {
                    anchors[1] = Some(q_pos + (targets[1] - ref_pos).max(0));
                }
                ref_pos += l;
                q_pos += l;
            }
            Cigar::Ins(l) => {
                q_pos += i64::from(*l);
            }
            Cigar::Del(l) | Cigar::RefSkip(l) => {
                let l = i64::from(*l);
                if anchors[0].is_none() && ref_pos + l >= targets[0] {
                    anchors[0] = Some(q_pos);
                }
                if anchors[1].is_none() && ref_pos + l > targets[1] {
                    anchors[1] = Some(q_pos);
                }
                ref_pos += l;
            }
            Cigar::SoftClip(l) => {
                // a clip reaching into the window means the read stops inside the locus
                if ref_pos >= targets[0] && ref_pos <= targets[1] {
                    return None;
                }
                q_pos += i64::from(*l);
            }
            Cigar::HardClip(_) => {
                if ref_pos >= targets[0] && ref_pos <= targets[1] {
                    return None;
                }
            }
            Cigar::Pad(_) => (),
        }
        if anchors[1].is_some() {
            break;
        }
    }
    // a read whose alignment ends exactly at the window edge still spans it
    if anchors[1].is_none() && ref_pos >= targets[1] {
        anchors[1] = Some(q_pos);
    }

    let (win_lo, win_hi) = (anchors[0]?, anchors[1]?);
    // trim as many query bases as the flank has reference bases: an insertion the aligner
    // parked in a flank lengthens the allele, a deletion there shortens it, and a flank
    // that matches the reference exactly leaves the annotated interval untouched
    let lo = win_lo + flank_lo;
    let hi = win_hi - flank_hi;
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
        // skip reads with mapq 0 or reads that do not span the repeat locus
        if r.mapq() == 0
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
    let _reads = get_overlapping_reads(&mut bam, &repeat, unphased, 60);
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
    let reads = get_overlapping_reads(&mut bam, &repeat, unphased, 60)
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
    fn flank_captures_an_insertion_shifted_out_of_the_interval() {
        // 5bp insertion at 0-based 996, four bases before the annotated repeat start
        let mut seq = vec![b'A'; 105];
        seq[46..51].copy_from_slice(b"TTTTT"); // the insertion
        seq[55..65].copy_from_slice(b"CGCGCGCGCG"); // 0-based 1000..1010
        let rec = record(950, vec![Cigar::Match(46), Cigar::Ins(5), Cigar::Match(54)], &seq);
        // without flank the insertion is missed and only the annotated interval is taken
        let tight = repeat_sequence_from_alignment(&rec, START, END, 0).unwrap();
        assert_eq!(tight, b"CGCGCGCGCG".to_vec());
        // with flank the allele grows by exactly the 5 inserted bases
        let padded = repeat_sequence_from_alignment(&rec, START, END, 10).unwrap();
        assert_eq!(padded.len(), 15);
        assert!(padded.ends_with(b"CGCGCGCGCG"));
    }

    #[test]
    fn flank_captures_a_deletion_shifted_out_of_the_interval() {
        // 4bp deletion at 0-based 996, just before the annotated repeat start
        let mut seq = vec![b'A'; 96];
        seq[46..56].copy_from_slice(b"CGCGCGCGCG"); // 0-based 1000..1010
        let rec = record(950, vec![Cigar::Match(46), Cigar::Del(4), Cigar::Match(50)], &seq);
        let tight = repeat_sequence_from_alignment(&rec, START, END, 0).unwrap();
        assert_eq!(tight, b"CGCGCGCGCG".to_vec());
        // with flank the allele shrinks by exactly the 4 deleted bases
        let padded = repeat_sequence_from_alignment(&rec, START, END, 10).unwrap();
        assert_eq!(padded.len(), 6);
    }

    #[test]
    fn flank_leaves_a_clean_locus_untouched() {
        let mut seq = vec![b'A'; 100];
        seq[50..60].copy_from_slice(b"CGCGCGCGCG");
        let rec = record(950, vec![Cigar::Match(100)], &seq);
        for flank in [0, 5, 10, 25] {
            let out = repeat_sequence_from_alignment(&rec, START, END, flank).unwrap();
            assert_eq!(out, b"CGCGCGCGCG".to_vec(), "flank {flank}");
        }
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
