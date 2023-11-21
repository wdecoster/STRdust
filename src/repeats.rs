use log::error;
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

    pub fn make_repeat_compressed_sequence(&self, fasta: &String, flanking: u32) -> Vec<u8> {
        let fas = faidx::Reader::from_path(fasta).expect("Failed to read fasta");
        let fai = format!("{}.fai", fasta);
        // check if the chromosome exists in the fai file
        if !std::fs::read_to_string(fai)
            .expect("Failed to read fai file")
            .lines()
            .any(|line| line.split('\t').collect::<Vec<&str>>()[0] == self.chrom)
        {
            error!("Chromosome {} not found in the fasta file", self.chrom);
            panic!();
        }
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
        // // write the new reference sequence to a file
        // let mut newref_file =
        //     std::fs::File::create("newref.fa").expect("Unable to create newref.fa");
        // newref_file
        //     .write_all(&[fas_left, fas_right].concat())
        //     .expect("Unable to write to newref.fa");
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

/// parse a region string
pub fn process_region(reg: &str) -> Result<RepeatInterval, Box<dyn std::error::Error>> {
    let chrom = reg.split(':').collect::<Vec<&str>>()[0].to_string();
    let interval = reg.split(':').collect::<Vec<&str>>()[1];
    let start: u32 = interval.split('-').collect::<Vec<&str>>()[0]
        .parse()
        .unwrap();
    let end: u32 = interval.split('-').collect::<Vec<&str>>()[1]
        .parse()
        .unwrap();
    assert!(
        end - start > 0,
        r#"Invalid region: begin has to be smaller than end."#
    );
    Ok(RepeatInterval { chrom, start, end })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_make_repeat_compressed_sequence() {
        let fasta = String::from("test_data/chr7.fa.gz");
        let repeat = RepeatInterval {
            chrom: String::from("chr7"),
            start: 154654404,
            end: 154654432,
        };
        let flanking = 2000;
        let newref = repeat.make_repeat_compressed_sequence(&fasta, flanking);
        // println!("{:?}", std::str::from_utf8(&seq[9990..10010]));
        // println!("{:?}", std::str::from_utf8(&seq));
        assert_eq!(newref.len() as u32, flanking * 2);
    }

    // captures the case when the repeat interval is out of bounds for the fasta file
    #[test]
    fn test_interval_out_of_bounds() {
        let fasta = String::from("test_data/chr7.fa.gz");
        let repeat = RepeatInterval {
            chrom: String::from("chr7"),
            start: 1154654404,
            end: 1154654432,
        };
        assert!(repeat.reference_repeat_sequence(&fasta).is_none());
    }

    // capture the case where the newref extended by <flanking> bases is out of bounds (<0) for the fasta file
    // and thus a shorter sequence is returned
    #[test]
    fn test_make_newref_interval_out_of_bounds1() {
        let fasta = String::from("test_data/chr7.fa.gz");
        let repeat = RepeatInterval {
            chrom: String::from("chr7"),
            start: 1000,
            end: 5010,
        };
        let flanking = 2000;
        let newref = repeat.make_repeat_compressed_sequence(&fasta, flanking);
        assert!((newref.len() as u32) < flanking * 2);
    }

    // capture the case where the newref extended by <flanking>b is out of bounds (> chromosome end) for the fasta file
    #[test]
    fn test_make_newref_interval_out_of_bounds2() {
        let fasta = String::from("test_data/chr7.fa.gz");
        let repeat = RepeatInterval {
            chrom: String::from("chr7"),
            start: 159345373,
            end: 159345400,
        };
        let flanking = 2000;
        let newref = repeat.make_repeat_compressed_sequence(&fasta, flanking);
        // println!("{:?}", std::str::from_utf8(&newref));
        // println!("{}", newref.len());
        assert!((newref.len() as u32) < flanking * 2);
    }

    #[test]
    #[should_panic]
    fn test_make_newref_interval_invalid_chrom() {
        let fasta = String::from("test_data/chr7.fa.gz");
        let repeat = RepeatInterval {
            chrom: String::from("chr27"),
            start: 159345373,
            end: 159345400,
        };
        let flanking = 2000;
        let _ = repeat.make_repeat_compressed_sequence(&fasta, flanking);
    }
}
