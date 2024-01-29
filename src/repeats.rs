use bio::io::bed;
use rust_htslib::faidx;
use std::fmt;

pub struct RepeatInterval {
    pub chrom: String,
    pub start: u32,
    pub end: u32,
}

impl fmt::Display for RepeatInterval {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}:{}-{}", self.chrom, self.start, self.end)
    }
}

impl RepeatInterval {
    /// Make a repeat compressed sequence from a fasta file
    /// The repeat sequence is removed, and <flanking> bases of flanking sequence is added on either side
    /// If the extended interval (+/- flanking) is out of bounds, the extended interval is shortened automatically by faidx

    // parse a bed record
    pub fn from_bed(rec: &bed::Record, fasta: &str) -> Option<Self> {
        let chrom = rec.chrom().to_string();
        let start = rec.start().try_into().unwrap();
        let end = rec.end().try_into().unwrap();
        RepeatInterval::new_interval(chrom, start, end, fasta)
    }

    // parse a region string
    pub fn from_string(reg: &str, fasta: &str) -> Option<Self> {
        let chrom = reg.split(':').collect::<Vec<&str>>()[0].to_string();
        let interval = reg.split(':').collect::<Vec<&str>>()[1];
        let start: u32 = interval.split('-').collect::<Vec<&str>>()[0]
            .parse()
            .unwrap();
        let end: u32 = interval.split('-').collect::<Vec<&str>>()[1]
            .parse()
            .unwrap();
        Self::new_interval(chrom, start, end, fasta)
    }

    fn new_interval(chrom: String, start: u32, end: u32, fasta: &str) -> Option<Self> {
        if end < start {
            return None;
        }

        let fai = format!("{}.fai", fasta);
        // check if the chromosome exists in the fai file
        // and if the end coordinate is within the chromosome length
        for line in std::fs::read_to_string(fai)
            .expect("Failed to read fai file")
            .lines()
        {
            if line.split('\t').collect::<Vec<&str>>()[0] == chrom
                && line.split('\t').collect::<Vec<&str>>()[1]
                    .parse::<u32>()
                    .expect("Failed parsing chromosome length from fai file")
                    > end
            {
                return Some(Self { chrom, start, end });
            }
        }
        // if the chromosome is not in the fai file or the end does not fit the interval, return None
        None
    }
    pub fn new(chrom: &str, start: u32, end: u32) -> Self {
        Self {
            chrom: chrom.to_string(),
            start,
            end,
        }
    }

    pub fn make_repeat_compressed_sequence(&self, fasta: &String, flanking: u32) -> Vec<u8> {
        let fas = faidx::Reader::from_path(fasta).expect("Failed to read fasta");
        let fas_left = fas
            .fetch_seq(
                &self.chrom,
                (self.start.saturating_sub(flanking + 2)) as usize,
                self.start as usize - 2,
            )
            .expect("Failed to extract fas_left sequence from fasta for {chrom}:{start}-{end}");
        let fas_right = fas
            .fetch_seq(
                &self.chrom,
                self.end as usize,
                (self.end + flanking - 2) as usize,
            )
            .expect("Failed to extract fas_right sequence from fasta for {chrom}:{start}-{end}");

        let newref = [fas_left, fas_right].concat();
        unsafe { libc::free(fas_left.as_ptr() as *mut std::ffi::c_void) }; // Free up memory (https://github.com/rust-bio/rust-htslib/issues/401#issuecomment-1704290171)
        unsafe { libc::free(fas_right.as_ptr() as *mut std::ffi::c_void) }; // Free up memory
        newref
    }

    pub fn reference_repeat_sequence(&self, fasta: &String) -> Option<String> {
        let fas = faidx::Reader::from_path(fasta).expect("Failed to read fasta");
        let repeat_ref_sequence = std::str::from_utf8(
            fas.fetch_seq(&self.chrom, self.start as usize - 1, self.end as usize)
                .expect("Failed to extract repeat sequence from fasta for {chrom}:{start}-{end}"),
        )
        .expect("Failed to convert repeat sequence to string for {chrom}:{start}-{end}")
        .to_string();

        // If the repeat sequence is out of bounds, None is returned
        if repeat_ref_sequence == "N" {
            eprintln!(
                "Cannot genotype repeat at {self} because it is out of bounds for the fasta file",
            );
            return None;
        }
        Some(repeat_ref_sequence)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_make_repeat_compressed_sequence() {
        let fasta = String::from("test_data/chr7.fa.gz");
        let repeat = RepeatInterval::from_string("chr7:154654404-154654432", &fasta).unwrap();
        let flanking = 2000;
        let newref = repeat.make_repeat_compressed_sequence(&fasta, flanking);
        assert_eq!(newref.len() as u32, flanking * 2);
    }

    // captures the case when the repeat interval is out of bounds for the fasta file
    #[test]
    fn test_interval_out_of_bounds() {
        let fasta = String::from("test_data/chr7.fa.gz");
        let repeat = RepeatInterval::from_string("chr7:1154654404-1154654432", &fasta);
        assert!(repeat.is_none());
    }

    // capture the case where the newref extended by <flanking> bases is out of bounds (<0) for the fasta file
    // and thus a shorter sequence is returned
    #[test]
    fn test_make_newref_interval_out_of_bounds1() {
        let fasta = String::from("test_data/chr7.fa.gz");
        let repeat = RepeatInterval::from_string("chr7:1000-5010", &fasta).unwrap();
        let flanking = 2000;
        let newref = repeat.make_repeat_compressed_sequence(&fasta, flanking);
        assert!((newref.len() as u32) < flanking * 2);
    }

    // capture the case where the newref extended by <flanking>b is out of bounds (> chromosome end) for the fasta file
    #[test]
    fn test_make_newref_interval_out_of_bounds2() {
        let fasta = String::from("test_data/chr7.fa.gz");
        let repeat = RepeatInterval::from_string("chr7:159345373-159345400", &fasta).unwrap();
        let flanking = 2000;
        let newref = repeat.make_repeat_compressed_sequence(&fasta, flanking);
        assert!((newref.len() as u32) < flanking * 2);
    }

    #[test]
    fn test_make_newref_interval_invalid_chrom() {
        let fasta = String::from("test_data/chr7.fa.gz");
        let repeat =
            crate::repeats::RepeatInterval::from_string("chr27:159345373-159345400", &fasta);
        assert!(repeat.is_none());
    }
}
