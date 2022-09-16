use flate2::read;
use std::ffi::OsStr;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Read normal or compressed files seamlessly
/// Uses the presence of a `.gz` extension to decide
pub fn reader(filename: &str) -> Box<dyn BufRead> {
    let path = Path::new(filename);
    let file = match File::open(&path) {
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

/// parse a region string
pub fn process_region(reg: String) -> Result<(String, u32, u32), Box<dyn std::error::Error>> {
    let chrom = reg.split(':').collect::<Vec<&str>>()[0];
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
    Ok((chrom.to_string(), start, end))
}
