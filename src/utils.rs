use flate2::read;
use log::error;
use std::ffi::OsStr;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

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

pub fn check_files_exist(args: &crate::Cli) {
    if !Path::new(&args.bam).exists() {
        error!("Alignment file not found: {}", args.bam);
        std::process::exit(1);
    }
    let index_extension = if args.bam.ends_with(".cram") {
        "crai"
    } else {
        "bai"
    };
    let index = format!("{}.{}", args.bam, index_extension);
    if !Path::new(&index).exists() {
        error!("Index file not found: {}", index);
        std::process::exit(1);
    }
    if !Path::new(&args.fasta).exists() {
        error!("FASTA file not found: {}", args.fasta);
        std::process::exit(1);
    }
    let fai = format!("{}.fai", args.fasta);
    if !Path::new(&fai).exists() {
        error!("FASTA index file not found: {}", fai);
        std::process::exit(1);
    }
}
