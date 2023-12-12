use rust_htslib::bam;
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::Read;
use std::collections::HashMap;
use url::Url;

pub struct Reads {
    // could consider not to use a hashmap here and use an attribute per phase
    pub seqs: HashMap<u8, Vec<Vec<u8>>>,
    pub ps: Option<u32>,
}

pub fn get_overlapping_reads(
    bamf: &String,
    repeat: &crate::repeats::RepeatInterval,
    unphased: bool,
) -> Option<Reads> {
    let mut bam = if bamf.starts_with("s3") || bamf.starts_with("https://") {
        bam::IndexedReader::from_url(&Url::parse(bamf.as_str()).expect("Failed to parse s3 URL"))
            .unwrap_or_else(|err| panic!("Error opening remote BAM: {err}"))
    } else {
        bam::IndexedReader::from_path(bamf)
            .unwrap_or_else(|err| panic!("Error opening local BAM: {err}"))
    };
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
        if r.mapq() < 10 {
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
    let repeat = crate::repeats::RepeatInterval {
        chrom: String::from("chr7"),
        start: 154654404,
        end: 154654432,
    };
    let unphased = false;
    let _reads = get_overlapping_reads(&bam, &repeat, unphased);
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
