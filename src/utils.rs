use flate2::read;
use rust_htslib::bam;
use rust_htslib::bam::record::Aux;
use std::ffi::OsStr;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

pub fn get_phase(record: &bam::Record) -> u8 {
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

/// Read normal or compressed files seamlessly
/// Uses the presence of a `.gz` extension to decide
pub fn reader(filename: &str) -> Box<dyn BufRead> {
    let path = Path::new(filename);
    let file = match File::open(path) {
        Err(why) => panic!("couldn't open {}: {}", path.display(), why),
        Ok(file) => file,
    };

    if path.extension() == Some(OsStr::new("gz")) {
        Box::new(BufReader::with_capacity(
            128 * 1024,
            read::GzDecoder::new(file),
        ))
    } else {
        Box::new(BufReader::with_capacity(128 * 1024, file))
    }
}

#[cfg(test)]
mod tests {
    use rust_htslib::bam::Read;

    use super::*;

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
}
