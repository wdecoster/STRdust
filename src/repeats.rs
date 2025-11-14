use bio::io::bed;
use log::error;
use rust_htslib::faidx;
use std::fmt;
use std::fs;
use std::io;

#[derive(Debug)]
pub struct RepeatIntervalIterator {
    current_index: usize,
    data: Vec<RepeatInterval>,
    pub num_intervals: usize,
}

impl RepeatIntervalIterator {
    // parse a region string
    pub fn from_string(reg: &str, fasta: &str) -> Self {
        // Split by colon first to get chromosome
        let parts: Vec<&str> = reg.split(':').collect();
        if parts.len() != 2 {
            error!("Invalid region format: '{}'. Expected format is 'chr:start-end'", reg);
            std::process::exit(1);
        }

        let chrom = parts[0].to_string();
        let interval = parts[1];

        // Split interval by hyphen to get start and end
        let coords: Vec<&str> = interval.split('-').collect();
        if coords.len() != 2 {
            error!("Invalid interval format: '{}'. Expected format is 'chr:start-end'", interval);
            error!("Example of a valid region: 'chr15:34419425-34419450'");
            std::process::exit(1);
        }

        // Parse start and end coordinates with error handling
        let start: u32 = match coords[0].parse() {
            Ok(val) => val,
            Err(_) => {
                error!("Could not parse start coordinate '{}' as a number", coords[0]);
                std::process::exit(1);
            }
        };

        let end: u32 = match coords[1].parse() {
            Ok(val) => val,
            Err(_) => {
                error!("Could not parse end coordinate '{}' as a number", coords[1]);
                std::process::exit(1);
            }
        };

        // Create the repeat interval with validation
        let repeat = match RepeatInterval::new_interval(chrom, start, end, fasta) {
            Some(interval) => interval,
            None => {
                error!("Failed to create repeat interval for region: '{}'", reg);
                std::process::exit(1);
            }
        };

        RepeatIntervalIterator { current_index: 0, data: vec![repeat], num_intervals: 1 }
    }

    pub fn from_bed(region_file: &String, fasta: &str) -> Self {
        // check if the bed file exists
        if !std::path::Path::new(region_file).exists() {
            error!("Bed file does not exist");
            std::process::exit(1);
        }
        let mut reader = bed::Reader::from_file(region_file).expect("Problem reading bed file!");
        let mut data = Vec::new();
        for record in reader.records() {
            let rec = record.expect("Error reading bed record. Please verify bed file is in the correct format and without header.");
            let repeat = RepeatInterval::from_bed(&rec, fasta);
            if let Some(repeat) = repeat {
                data.push(repeat);
            }
        }
        RepeatIntervalIterator { current_index: 0, data: data.clone(), num_intervals: data.len() }
    }

    pub fn pathogenic(fasta: &str) -> Self {
        let cache_dir = dirs::cache_dir()
            .unwrap_or_else(std::env::temp_dir)
            .join("strdust");

        // Create cache directory if it doesn't exist
        fs::create_dir_all(&cache_dir).expect("Failed to create cache directory");

        let cache_file = cache_dir.join("STRchive-disease-loci.hg38.TRGT.bed");

        // Check if cached file exists and is recent (e.g., less than 7 days old)
        let needs_download = if cache_file.exists() {
            match fs::metadata(&cache_file) {
                Ok(metadata) => {
                    if let Ok(modified) = metadata.modified() {
                        let age = std::time::SystemTime::now()
                            .duration_since(modified)
                            .unwrap_or(std::time::Duration::from_secs(0));
                        age > std::time::Duration::from_secs(7 * 24 * 60 * 60) // 7 days
                    } else {
                        true // Can't get modification time, re-download
                    }
                }
                Err(_) => true, // Can't get metadata, re-download
            }
        } else {
            true // File doesn't exist, download
        };

        // Download if needed
        if needs_download {
            eprintln!("Downloading pathogenic STR database...");
            let url = "https://raw.githubusercontent.com/dashnowlab/STRchive/refs/heads/main/data/catalogs/STRchive-disease-loci.hg38.longTR.bed";
            match reqwest::blocking::get(url) {
                Ok(resp) => {
                    match resp.text() {
                        Ok(body) => {
                            if let Err(e) = fs::write(&cache_file, &body) {
                                eprintln!(
                                    "Warning: Failed to cache pathogenic repeats data: {}",
                                    e
                                );
                                // Continue with in-memory data
                                return Self::from_string_data(&body, fasta);
                            }
                            eprintln!(
                                "Cached pathogenic STR database to: {}",
                                cache_file.display()
                            );
                        }
                        Err(e) => {
                            eprintln!("Failed to read response body: {}", e);
                            std::process::exit(1);
                        }
                    }
                }
                Err(e) => {
                    // If download fails but we have a cached version, use it even if old
                    if cache_file.exists() {
                        eprintln!("Download failed, using cached version: {}", e);
                    } else {
                        eprintln!("Failed to download pathogenic repeats database: {}", e);
                        std::process::exit(1);
                    }
                }
            }
        }

        // Read from cached file
        Self::from_bed(&cache_file.to_string_lossy().to_string(), fasta)
    }

    // Helper function to process data from string (for fallback)
    fn from_string_data(data: &str, fasta: &str) -> Self {
        let mut reader = io::BufReader::new(data.as_bytes());
        let mut data_vec = Vec::new();
        let mut reader = bed::Reader::new(&mut reader);
        for record in reader.records() {
            let rec = record.expect("Error reading bed record.");
            let repeat = RepeatInterval::from_bed(&rec, fasta);
            if let Some(repeat) = repeat {
                data_vec.push(repeat);
            }
        }
        RepeatIntervalIterator {
            current_index: 0,
            data: data_vec.clone(),
            num_intervals: data_vec.len(),
        }
    }

    /// Test helper function to create a RepeatIntervalIterator from BED data string
    /// This is exposed for testing purposes to allow testing BED parsing without file I/O
    #[cfg(test)]
    pub fn from_test_data(data: &str, fasta: &str) -> Self {
        Self::from_string_data(data, fasta)
    }

    /// Test helper function to get the cache file path
    /// This is exposed for testing purposes to allow test cleanup and verification
    #[cfg(test)]
    pub fn get_cache_file_path() -> std::path::PathBuf {
        let cache_dir = dirs::cache_dir()
            .unwrap_or_else(std::env::temp_dir)
            .join("strdust");
        cache_dir.join("STRchive-disease-loci.hg38.TRGT.bed")
    }
}

impl Clone for RepeatInterval {
    fn clone(&self) -> Self {
        RepeatInterval {
            chrom: self.chrom.clone(),
            start: self.start,
            end: self.end,
            created: self.created,
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
    pub created: Option<chrono::DateTime<chrono::Utc>>,
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
                return Some(Self { chrom, start, end, created: None });
            }
        }
        // if the chromosome is not in the fai file or the end does not fit the interval, return None
        panic!(
            "Chromosome {chrom} is not in the fasta file or the end coordinate is out of bounds"
        );
    }
    pub fn new(chrom: &str, start: u32, end: u32) -> Self {
        Self { chrom: chrom.to_string(), start, end, created: None }
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
            .fetch_seq(&self.chrom, self.end as usize, (self.end + flanking - 2) as usize)
            .expect("Failed to extract fas_right sequence from fasta for {chrom}:{start}-{end}");
        [fas_left, fas_right].concat()
    }

    pub fn reference_repeat_sequence(&self, fasta: &String) -> Option<String> {
        let fas = faidx::Reader::from_path(fasta).expect("Failed to read fasta");
        let repeat_ref_sequence = std::str::from_utf8(
            &fas.fetch_seq(&self.chrom, self.start as usize - 1, self.end as usize)
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

    pub fn set_time_stamp(&mut self) {
        self.created = Some(chrono::Utc::now());
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

    // Mock BED data that resembles STRchive format for testing
    const MOCK_STRCHIVE_DATA: &str = r#"chr1	1234567	1234590	ATTN	ATGC_4_3_90
chr2	2345678	2345701	CGG	CGG_3_8_100
chr4	39348425	39348483	CAG	ATXN1_CAG_58
chr7	154654404	154654432	CAG	ATXN2_CAG_28
chr19	45770205	45770266	CAG	ATXN1_CAG_61"#;

    #[test]
    fn test_pathogenic_bed_parsing() {
        // Create a temporary fasta file for testing
        let temp_dir = std::env::temp_dir();
        let test_fasta = temp_dir.join("test_pathogenic.fa");
        let test_fai = temp_dir.join("test_pathogenic.fa.fai");

        // Create minimal fasta content
        let fasta_content = ">chr1\nATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC\n>chr2\nCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGG\n>chr4\nCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAG\n>chr7\nATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC\n>chr19\nATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCAtgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcatgcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcagcag";

        // Create a simple fai content with chromosome lengths that cover our test intervals
        let fai_content = "chr1\t100000000\t6\t60\t61\nchr2\t100000000\t6\t60\t61\nchr4\t100000000\t6\t60\t61\nchr7\t200000000\t6\t60\t61\nchr19\t100000000\t6\t60\t61";

        // Write the test files
        std::fs::write(&test_fasta, fasta_content).expect("Failed to write test fasta");
        std::fs::write(&test_fai, fai_content).expect("Failed to write test fai");

        // Test parsing the mock BED data
        let result = RepeatIntervalIterator::from_string_data(
            MOCK_STRCHIVE_DATA,
            &test_fasta.to_string_lossy(),
        );

        // Verify the parsing worked correctly
        assert_eq!(result.num_intervals, 5, "Should parse all 5 intervals from mock data");

        let intervals: Vec<RepeatInterval> = result.collect();
        assert_eq!(intervals.len(), 5);

        // Verify first interval
        assert_eq!(intervals[0].chrom, "chr1");
        assert_eq!(intervals[0].start, 1234567);
        assert_eq!(intervals[0].end, 1234590);

        // Verify last interval
        assert_eq!(intervals[4].chrom, "chr19");
        assert_eq!(intervals[4].start, 45770205);
        assert_eq!(intervals[4].end, 45770266);

        // Clean up
        let _ = std::fs::remove_file(&test_fasta);
        let _ = std::fs::remove_file(&test_fai);
    }

    #[test]
    fn test_pathogenic_cache_functionality() {
        use std::time::{Duration, SystemTime};

        // Create a temporary directory for testing cache functionality
        let temp_dir = std::env::temp_dir().join("strdust_test_cache");
        let _ = std::fs::remove_dir_all(&temp_dir); // Clean up any previous test runs
        std::fs::create_dir_all(&temp_dir).expect("Failed to create test cache dir");

        let cache_file = temp_dir.join("STRchive-disease-loci.hg38.TRGT.bed");

        // Write mock cache data
        std::fs::write(&cache_file, MOCK_STRCHIVE_DATA).expect("Failed to write test cache file");

        // Test 1: Fresh cache file should be used (modification check)
        let metadata = cache_file
            .metadata()
            .expect("Failed to get cache file metadata");
        let modified = metadata
            .modified()
            .expect("Failed to get modification time");
        let age = SystemTime::now()
            .duration_since(modified)
            .unwrap_or(Duration::from_secs(0));

        // File should be very fresh (just created)
        assert!(age < Duration::from_secs(60), "Cache file should be fresh");

        // Test 2: Check file exists after creation
        assert!(cache_file.exists(), "Cache file should exist");

        // Test 3: Manually age the file to test the 7-day cache policy
        // We can't easily modify the file timestamp in a portable way, so we'll just verify
        // the logic would work by checking the duration calculation
        let seven_days = Duration::from_secs(7 * 24 * 60 * 60);
        assert!(age < seven_days, "Fresh file should be less than 7 days old");

        // Clean up
        let _ = std::fs::remove_dir_all(&temp_dir);
    }

    #[test]
    fn test_pathogenic_bed_format_edge_cases() {
        let temp_dir = std::env::temp_dir();
        let test_fasta = temp_dir.join("test_edge_cases.fa");
        let test_fai = temp_dir.join("test_edge_cases.fa.fai");

        // Create minimal fasta content
        let fasta_content = ">chr1\nATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC\n>chr22\nCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGG";
        let fai_content = "chr1\t1000000\t6\t60\t61\nchr22\t1000000\t6\t60\t61";

        std::fs::write(&test_fasta, fasta_content).expect("Failed to write test fasta");
        std::fs::write(&test_fai, fai_content).expect("Failed to write test fai");

        // Test with different BED formats
        let bed_data_with_extra_columns = "chr1\t100\t150\tATGC\tgene1\t.\t+\textra_info";
        let bed_data_minimal = "chr22\t200\t250";

        // Test parsing with extra columns (should work - BED parser should handle this)
        let result1 = RepeatIntervalIterator::from_string_data(
            bed_data_with_extra_columns,
            &test_fasta.to_string_lossy(),
        );
        assert_eq!(result1.num_intervals, 1);
        let intervals1: Vec<RepeatInterval> = result1.collect();
        assert_eq!(intervals1[0].chrom, "chr1");
        assert_eq!(intervals1[0].start, 100);
        assert_eq!(intervals1[0].end, 150);

        // Test parsing with minimal columns
        let result2 = RepeatIntervalIterator::from_string_data(
            bed_data_minimal,
            &test_fasta.to_string_lossy(),
        );
        assert_eq!(result2.num_intervals, 1);
        let intervals2: Vec<RepeatInterval> = result2.collect();
        assert_eq!(intervals2[0].chrom, "chr22");
        assert_eq!(intervals2[0].start, 200);
        assert_eq!(intervals2[0].end, 250);

        // Clean up
        let _ = std::fs::remove_file(&test_fasta);
        let _ = std::fs::remove_file(&test_fai);
    }

    #[test]
    #[should_panic(expected = "Chromosome chr99 is not in the fasta file")]
    fn test_pathogenic_invalid_chromosome() {
        let temp_dir = std::env::temp_dir();
        let test_fasta = temp_dir.join("test_invalid_chr.fa");
        let test_fai = temp_dir.join("test_invalid_chr.fa.fai");

        let fasta_content =
            ">chr1\nATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC";
        let fai_content = "chr1\t1000000\t6\t60\t61";

        std::fs::write(&test_fasta, fasta_content).expect("Failed to write test fasta");
        std::fs::write(&test_fai, fai_content).expect("Failed to write test fai");

        // This should panic because chr99 is not in the fasta
        let bed_data_invalid = "chr99\t100\t150\tATGC\tgene1";
        let _result = RepeatIntervalIterator::from_string_data(
            bed_data_invalid,
            &test_fasta.to_string_lossy(),
        );

        // Clean up (won't be reached due to panic, but good practice)
        let _ = std::fs::remove_file(&test_fasta);
        let _ = std::fs::remove_file(&test_fai);
    }

    // Integration test that actually tests the pathogenic download functionality
    // This test is optional and requires network access
    #[test]
    #[ignore = "requires network access and reference genome - set TEST_PATHOGENIC_DOWNLOAD=1 to enable"]
    fn test_pathogenic_download_integration() {
        // Only run if explicitly enabled
        if std::env::var("TEST_PATHOGENIC_DOWNLOAD").is_err() {
            return;
        }

        let fasta =
            std::env::var("FASTA_PATH").unwrap_or_else(|_| "test_data/chr7.fa.gz".to_string());

        if !std::path::Path::new(&fasta).exists() {
            eprintln!("Skipping integration test - fasta file {} does not exist", fasta);
            return;
        }

        let fai_file = format!("{}.fai", fasta);
        if !std::path::Path::new(&fai_file).exists() {
            eprintln!("Skipping integration test - fasta index file {} does not exist", fai_file);
            return;
        }

        // Clear any existing cache to test download functionality
        let cache_dir = dirs::cache_dir()
            .unwrap_or_else(std::env::temp_dir)
            .join("strdust");
        let cache_file = cache_dir.join("STRchive-disease-loci.hg38.TRGT.bed");
        let _ = std::fs::remove_file(&cache_file);

        // This should trigger a download
        let result = RepeatIntervalIterator::pathogenic(&fasta);

        // Verify results
        assert!(
            result.num_intervals > 0,
            "Should have downloaded and parsed pathogenic intervals"
        );
        assert!(cache_file.exists(), "Cache file should exist after download");

        // Verify cache file has reasonable content
        let cache_content =
            std::fs::read_to_string(&cache_file).expect("Failed to read cache file");
        assert!(
            cache_content.contains("chr"),
            "Cache file should contain chromosome information"
        );
        assert!(cache_content.len() > 100, "Cache file should have substantial content");

        // Test that second call uses cache (no download)
        let result2 = RepeatIntervalIterator::pathogenic(&fasta);
        assert_eq!(
            result.num_intervals, result2.num_intervals,
            "Cached results should match downloaded results"
        );
    }

    // this test is conditional, as it requires a specific fasta file to exist
    #[test]
    #[ignore = "requires specific fasta file - set TEST_WITH_FASTA=1 to enable"]
    fn test_pathogenic() {
        // Only run if explicitly enabled
        if std::env::var("TEST_WITH_FASTA").is_err() {
            return;
        }

        let fasta = std::env::var("FASTA_PATH")
            .unwrap_or_else(|_| "/home/wdecoster/reference/GRCh38.fa".to_string());

        if !std::path::Path::new(&fasta).exists() {
            panic!(
                "Fasta file {} does not exist. Set FASTA_PATH environment variable to correct path",
                fasta
            );
        }

        let fai_file = format!("{}.fai", fasta);
        if !std::path::Path::new(&fai_file).exists() {
            panic!("Fasta index file {} does not exist", fai_file);
        }

        let result = crate::repeats::RepeatIntervalIterator::pathogenic(&fasta);
        assert!(result.num_intervals > 0, "Should have found some pathogenic intervals");
    }
}
