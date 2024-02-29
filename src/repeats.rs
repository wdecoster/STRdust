use bio::io::bed;
use rust_htslib::faidx;
use std::fmt;
use std::io;

#[derive(Debug)]
pub struct RepeatIntervalIterator {
    // Add fields here that you need for iteration.
    // For example, you might need a current index and a vector of RepeatInterval.
    current_index: usize,
    data: Vec<RepeatInterval>,
}

impl RepeatIntervalIterator {
    // parse a region string
    pub fn from_string(reg: &str, fasta: &str) -> Self {
        let chrom = reg.split(':').collect::<Vec<&str>>()[0].to_string();
        let interval = reg.split(':').collect::<Vec<&str>>()[1];
        let start: u32 = interval.split('-').collect::<Vec<&str>>()[0]
            .parse()
            .unwrap();
        let end: u32 = interval.split('-').collect::<Vec<&str>>()[1]
            .parse()
            .unwrap();
        let repeat = RepeatInterval::new_interval(chrom, start, end, fasta)
            .expect("Failed to create repeat interval");
        RepeatIntervalIterator {
            current_index: 0,
            data: vec![repeat],
        }
    }
    pub fn from_bed(region_file: &String, fasta: &str) -> Self {
        let mut reader = bed::Reader::from_file(region_file).expect("Problem reading bed file!");
        let mut data = Vec::new();
        for record in reader.records() {
            let rec = record.expect("Error reading bed record.");
            let repeat = RepeatInterval::from_bed(&rec, fasta);
            if let Some(repeat) = repeat {
                data.push(repeat);
            }
        }
        RepeatIntervalIterator {
            current_index: 0,
            data,
        }
    }

    pub fn pathogenic(fasta: &str) -> Self {
        let url = "https://raw.githubusercontent.com/hdashnow/STRchive/main/data/hg38.STRchive-disease-loci.TRGT.bed";
        let resp = reqwest::blocking::get(url).expect("request to STRchive failed");
        let body = resp.text().expect("body invalid");
        let mut reader = io::BufReader::new(body.as_bytes());
        let mut data = Vec::new();
        let mut reader = bed::Reader::new(&mut reader);
        for record in reader.records() {
            let rec = record.expect("Error reading bed record.");
            let repeat = RepeatInterval::from_bed(&rec, fasta);
            if let Some(repeat) = repeat {
                data.push(repeat);
            }
        }
        RepeatIntervalIterator {
            current_index: 0,
            data,
        }
    }
}

impl Clone for RepeatInterval {
    fn clone(&self) -> Self {
        RepeatInterval {
            chrom: self.chrom.clone(),
            start: self.start,
            end: self.end,
        }
    }
}

impl Iterator for RepeatIntervalIterator {
    type Item = RepeatInterval;

    fn next(&mut self) -> Option<Self::Item> {
        // Implement the logic to get the next RepeatInterval here.
        // This is a simple example that gets the next item from a vector.
        if self.current_index < self.data.len() {
            let result = self.data[self.current_index].clone();
            self.current_index += 1;
            Some(result)
        } else {
            None
        }
    }
}

#[derive(Debug)]
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

    fn new_interval(chrom: String, start: u32, end: u32, fasta: &str) -> Option<Self> {
        if end < start {
            panic!("End coordinate is smaller than start coordinate for {chrom}:{start}-{end}")
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
        panic!(
            "Chromosome {chrom} is not in the fasta file or the end coordinate is out of bounds"
        );
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
        let repeat = RepeatIntervalIterator::from_string("chr7:154654404-154654432", &fasta);
        let flanking = 2000;
        for r in repeat {
            let newref = r.make_repeat_compressed_sequence(&fasta, flanking);
            assert_eq!(newref.len() as u32, flanking * 2);
        }
    }

    // captures the case when the repeat interval is out of bounds for the fasta file
    #[test]
    #[should_panic]
    fn test_interval_out_of_bounds() {
        let fasta = String::from("test_data/chr7.fa.gz");
        let _ = RepeatIntervalIterator::from_string("chr7:1154654404-1154654432", &fasta);
    }

    // capture the case where the newref extended by <flanking> bases is out of bounds (<0) for the fasta file
    // and thus a shorter sequence is returned
    #[test]
    fn test_make_newref_interval_out_of_bounds1() {
        let fasta = String::from("test_data/chr7.fa.gz");
        let repeat = RepeatIntervalIterator::from_string("chr7:1000-5010", &fasta);
        let flanking = 2000;
        for r in repeat {
            let newref = r.make_repeat_compressed_sequence(&fasta, flanking);
            assert!((newref.len() as u32) < flanking * 2);
        }
    }

    // capture the case where the newref extended by <flanking>b is out of bounds (> chromosome end) for the fasta file
    #[test]
    fn test_make_newref_interval_out_of_bounds2() {
        let fasta = String::from("test_data/chr7.fa.gz");
        let repeat = RepeatIntervalIterator::from_string("chr7:159345373-159345400", &fasta);
        let flanking = 2000;
        for r in repeat {
            let newref = r.make_repeat_compressed_sequence(&fasta, flanking);
            assert!((newref.len() as u32) < flanking * 2);
        }
    }

    #[test]
    #[should_panic]
    fn test_make_newref_interval_invalid_chrom() {
        let fasta = String::from("test_data/chr7.fa.gz");
        let _ = crate::repeats::RepeatIntervalIterator::from_string(
            "chr27:159345373-159345400",
            &fasta,
        );
    }

    // this test is ignored as it uses a file outside the test_data directory
    #[test]
    #[ignore]
    fn test_pathogenic() {
        let fasta = String::from("/home/wdecoster/reference/GRCh38.fa");
        let _ = crate::repeats::RepeatIntervalIterator::pathogenic(&fasta);
    }
}
