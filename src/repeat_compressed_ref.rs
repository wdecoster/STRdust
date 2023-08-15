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
    chrom: &String,
    start: u32,
    end: u32,
) -> Option<(Vec<u8>, String)> {
    let fas = faidx::Reader::from_path(fasta).expect("Failed to read fasta");
    let fai = format!("{}.fai", fasta);
    // check if the chromosome exists in the fai file
    if !std::fs::read_to_string(fai)
        .expect("Failed to read fai file")
        .lines()
        .any(|line| line.split('\t').collect::<Vec<&str>>()[0] == chrom)
    {
        error!("Chromosome {} not found in the fasta file", chrom);
        panic!();
    }
    let fas_left = fas
        .fetch_seq(
            chrom,
            (start.saturating_sub(10002)) as usize,
            start as usize - 2,
        )
        .expect("Failed to extract fas_left sequence from fasta for {chrom}:{start}-{end}");
    let fas_right = fas
        .fetch_seq(chrom, end as usize, (end + 9998) as usize)
        .expect("Failed to extract fas_right sequence from fasta for {chrom}:{start}-{end}");

    let repeat_ref_sequence = std::str::from_utf8(
        fas.fetch_seq(chrom, start as usize - 1, end as usize)
            .expect("Failed to extract repeat sequence from fasta for {chrom}:{start}-{end}"),
    )
    .expect("Failed to convert repeat sequence to string for {chrom}:{start}-{end}")
    .to_string();
    if repeat_ref_sequence == "N" {
        eprintln!(
            "Cannot genotype repeat at {chrom}:{start}-{end} because it is out of bounds for the fasta file",
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
        let chrom = String::from("chr7");
        let start = 154654404;
        let end = 154654432;
        let (newref, _) = make_repeat_compressed_sequence(&fasta, &chrom, start, end)
            .expect("Unable to make repeat compressed sequence");
        // println!("{:?}", std::str::from_utf8(&seq[9990..10010]));
        // println!("{:?}", std::str::from_utf8(&seq));
        assert_eq!(newref.len(), 20000);
    }

    // captures the case when the repeat interval is out of bounds for the fasta file
    #[test]
    fn test_interval_out_of_bounds() {
        let fasta = String::from("test_data/chr7.fa.gz");
        let chrom = String::from("chr7");
        let start = 1154654404;
        let end = 1154654432;
        assert!(make_repeat_compressed_sequence(&fasta, &chrom, start, end).is_none());
    }

    // capture the case where the newref extended by 10kb is out of bounds (<0) for the fasta file
    #[test]
    fn test_make_newref_interval_out_of_bounds1() {
        let fasta = String::from("test_data/chr7.fa.gz");
        let chrom = String::from("chr7");
        let start = 5000;
        let end = 5010;
        let (newref, _) = make_repeat_compressed_sequence(&fasta, &chrom, start, end)
            .expect("Unable to make repeat compressed sequence");
        assert!(newref.len() < 20000);
    }

    // capture the case where the newref extended by 10kb is out of bounds (> chromosome end) for the fasta file
    #[test]
    fn test_make_newref_interval_out_of_bounds2() {
        let fasta = String::from("test_data/chr7.fa.gz");
        let chrom = String::from("chr7");
        let start = 159345373;
        let end = 159345400;
        let (newref, _) = make_repeat_compressed_sequence(&fasta, &chrom, start, end)
            .expect("Unable to make repeat compressed sequence");
        // println!("{:?}", std::str::from_utf8(&newref));
        // println!("{}", newref.len());
        assert!(newref.len() < 20000);
    }

    #[test]
    #[should_panic]
    fn test_make_newref_interval_invalid_chrom() {
        let fasta = String::from("test_data/chr7.fa.gz");
        let chrom = String::from("chr27");
        let start = 159345373;
        let end = 159345400;
        let (_, _) = make_repeat_compressed_sequence(&fasta, &chrom, start, end)
            .expect("Unable to make repeat compressed sequence");
    }
}
