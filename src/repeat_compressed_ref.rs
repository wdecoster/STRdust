use std::io::Write;

use log::error;
use rust_htslib::faidx;

/// Make a repeat compressed sequence from a fasta file
/// The repeat sequence is removed, and 10kb of flanking sequence is added on either side
/// The repeat sequence is returned separately
/// If the repeat sequence is out of bounds, None is returned
/// If the extended interval (+/- 10kb) is out of bounds, the extended interval is shortened automatically by faidx
pub fn make_repeat_compressed_sequence(
    fasta: &String,
    repeat: &crate::utils::RepeatInterval,
    flanking: u32,
) -> Option<(Vec<u8>, String)> {
    let fas = faidx::Reader::from_path(fasta).expect("Failed to read fasta");
    let fai = format!("{}.fai", fasta);
    // check if the chromosome exists in the fai file
    if !std::fs::read_to_string(fai)
        .expect("Failed to read fai file")
        .lines()
        .any(|line| line.split('\t').collect::<Vec<&str>>()[0] == repeat.chrom)
    {
        error!("Chromosome {} not found in the fasta file", repeat.chrom);
        panic!();
    }
    let fas_left = fas
        .fetch_seq(
            &repeat.chrom,
            (repeat.start.saturating_sub(flanking + 2)) as usize,
            repeat.start as usize - 2,
        )
        .expect("Failed to extract fas_left sequence from fasta for {chrom}:{start}-{end}");
    let fas_right = fas
        .fetch_seq(
            &repeat.chrom,
            repeat.end as usize,
            (repeat.end + flanking - 2) as usize,
        )
        .expect("Failed to extract fas_right sequence from fasta for {chrom}:{start}-{end}");

    let repeat_ref_sequence = std::str::from_utf8(
        fas.fetch_seq(
            &repeat.chrom,
            repeat.start as usize - 1,
            repeat.end as usize,
        )
        .expect("Failed to extract repeat sequence from fasta for {chrom}:{start}-{end}"),
    )
    .expect("Failed to convert repeat sequence to string for {chrom}:{start}-{end}")
    .to_string();
    if repeat_ref_sequence == "N" {
        eprintln!(
            "Cannot genotype repeat at {repeat} because it is out of bounds for the fasta file",
        );
        return None;
    }
    // write the new reference sequence to a file
    let mut newref_file = std::fs::File::create("newref.fa").expect("Unable to create newref.fa");
    newref_file
        .write_all(&[fas_left, fas_right].concat())
        .expect("Unable to write to newref.fa");
    Some(([fas_left, fas_right].concat(), repeat_ref_sequence))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_make_repeat_compressed_sequence() {
        let fasta = String::from("test_data/chr7.fa.gz");
        let repeat = crate::utils::RepeatInterval {
            chrom: String::from("chr7"),
            start: 154654404,
            end: 154654432,
        };
        let flanking = 2000;
        let (newref, _) = make_repeat_compressed_sequence(&fasta, &repeat, flanking)
            .expect("Unable to make repeat compressed sequence");
        // println!("{:?}", std::str::from_utf8(&seq[9990..10010]));
        // println!("{:?}", std::str::from_utf8(&seq));
        assert_eq!(newref.len() as u32, flanking * 2);
    }

    // captures the case when the repeat interval is out of bounds for the fasta file
    #[test]
    fn test_interval_out_of_bounds() {
        let fasta = String::from("test_data/chr7.fa.gz");
        let repeat = crate::utils::RepeatInterval {
            chrom: String::from("chr7"),
            start: 1154654404,
            end: 1154654432,
        };
        let flanking = 2000;
        assert!(make_repeat_compressed_sequence(&fasta, &repeat, flanking).is_none());
    }

    // capture the case where the newref extended by 10kb is out of bounds (<0) for the fasta file
    #[test]
    fn test_make_newref_interval_out_of_bounds1() {
        let fasta = String::from("test_data/chr7.fa.gz");
        let repeat = crate::utils::RepeatInterval {
            chrom: String::from("chr7"),
            start: 1000,
            end: 5010,
        };
        let flanking = 2000;
        let (newref, _) = make_repeat_compressed_sequence(&fasta, &repeat, flanking)
            .expect("Unable to make repeat compressed sequence");
        assert!((newref.len() as u32) < flanking * 2);
    }

    // capture the case where the newref extended by 10kb is out of bounds (> chromosome end) for the fasta file
    #[test]
    fn test_make_newref_interval_out_of_bounds2() {
        let fasta = String::from("test_data/chr7.fa.gz");
        let repeat = crate::utils::RepeatInterval {
            chrom: String::from("chr7"),
            start: 159345373,
            end: 159345400,
        };
        let flanking = 2000;
        let (newref, _) = make_repeat_compressed_sequence(&fasta, &repeat, flanking)
            .expect("Unable to make repeat compressed sequence");
        // println!("{:?}", std::str::from_utf8(&newref));
        // println!("{}", newref.len());
        assert!((newref.len() as u32) < flanking * 2);
    }

    #[test]
    #[should_panic]
    fn test_make_newref_interval_invalid_chrom() {
        let fasta = String::from("test_data/chr7.fa.gz");
        let repeat = crate::utils::RepeatInterval {
            chrom: String::from("chr27"),
            start: 159345373,
            end: 159345400,
        };
        let flanking = 2000;
        let (_, _) = make_repeat_compressed_sequence(&fasta, &repeat, flanking)
            .expect("Unable to make repeat compressed sequence");
    }
}
