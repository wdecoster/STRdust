use log::debug;
use rust_htslib::bam;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::Read;
use std::collections::HashMap;
use std::env;
use url::Url;

pub struct Reads {
    // could consider not to use a hashmap here and use an attribute per phase
    pub seqs: HashMap<u8, Vec<Vec<u8>>>,
    pub ps: Option<u32>,
}

pub fn create_bam_reader(bamf: &str, fasta: &str) -> bam::IndexedReader {
    let mut bam = if bamf.starts_with("s3") || bamf.starts_with("https://") {
        if env::var("CURL_CA_BUNDLE").is_err() {
            env::set_var("CURL_CA_BUNDLE", "/etc/ssl/certs/ca-certificates.crt");
        }
        bam::IndexedReader::from_url(&Url::parse(bamf).expect("Failed to parse s3 URL"))
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

pub fn get_overlapping_reads(
    bam: &mut bam::IndexedReader,
    repeat: &crate::repeats::RepeatInterval,
    unphased: bool,
) -> Option<Reads> {
    let tid = bam
        .header()
        .tid(repeat.chrom.as_bytes())
        .unwrap_or_else(|| panic!("Invalid chromosome {}", repeat.chrom));
    bam.fetch((tid, repeat.start, repeat.end))
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
            debug!(
                "Skipping read {}",
                std::str::from_utf8(r.qname()).expect("Could get read identifier")
            );
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
            eprintln!("Cannot genotype {repeat}: no reads found");
        } else {
            eprintln!("Cannot genotype {repeat}: no phased reads found");
        }
        None
    } else {
        Some(Reads { seqs, ps })
    }
}

fn get_phase(record: &bam::Record) -> u8 {
    match record.aux(b"HP") {
        Ok(value) => {
            if let Aux::U8(v) = value {
                v
            } else {
                panic!("Unexpected type of Aux {value:?}")
            }
        }
        Err(_e) => 0,
    }
}

fn get_phase_set(record: &bam::Record) -> Option<u32> {
    match record.aux(b"PS") {
        Ok(value) => {
            if let Aux::U32(v) = value {
                Some(v)
            } else {
                panic!("Unexpected type of Aux {value:?}")
            }
        }
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
    };
    let unphased = false;
    let mut bam = create_bam_reader(&bam, &fasta);
    let _reads = get_overlapping_reads(&mut bam, &repeat, unphased);
}

#[test]
fn test_get_overlapping_reads_url1() {
    let bam = String::from("https://s3.amazonaws.com/1000g-ont/FIRST_100_FREEZE/minimap2_2.24_alignment_data/GM18501/GM18501.LSK110.R9.guppy646.sup.with5mC.pass.phased.bam");
    let fasta = String::from("test_data/chr7.fa.gz");
    let repeat = crate::repeats::RepeatInterval {
        chrom: String::from("chr20"),
        start: 154654404,
        end: 154654432,
    };
    let unphased = false;
    let mut bam = create_bam_reader(&bam, &fasta);
    let _reads = get_overlapping_reads(&mut bam, &repeat, unphased);
}

#[test]
fn test_get_overlapping_reads_url2() {
    let bam = String::from("s3://1000g-ont/FIRST_100_FREEZE/minimap2_2.24_alignment_data/GM18501/GM18501.LSK110.R9.guppy646.sup.with5mC.pass.phased.bam");
    let fasta = String::from("test_data/chr7.fa.gz");
    let repeat = crate::repeats::RepeatInterval {
        chrom: String::from("chr20"),
        start: 154654404,
        end: 154654432,
    };
    let unphased = false;
    let mut bam = create_bam_reader(&bam, &fasta);
    let _reads = get_overlapping_reads(&mut bam, &repeat, unphased);
}

#[test]
fn test_get_overlapping_reads_url3() {
    let bam = String::from("https://s3.amazonaws.com/1000g-ont/FIRST_100_FREEZE/minimap2_2.24_alignment_data/GM18501/GM18501.LSK110.R9.guppy646.sup.with5mC.pass.phased.bam");
    let fasta = String::from("test_data/chr7.fa.gz");
    let repeat = crate::repeats::RepeatInterval {
        chrom: String::from("chr20"),
        start: 154654404,
        end: 154654432,
    };
    let unphased = false;
    let mut bam = create_bam_reader(&bam, &fasta);
    let _reads = get_overlapping_reads(&mut bam, &repeat, unphased);
}

#[test]
fn test_get_overlapping_reads_url4() {
    let bam = String::from("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/hg38/HG00096.hg38.cram");
    let fasta = String::from("test_data/chr7.fa.gz");
    let repeat = crate::repeats::RepeatInterval {
        chrom: String::from("chr20"),
        start: 154654404,
        end: 154654432,
    };
    let unphased = false;
    let mut bam = create_bam_reader(&bam, &fasta);
    let _reads = get_overlapping_reads(&mut bam, &repeat, unphased);
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
