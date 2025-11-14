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
        Box::new(BufReader::with_capacity(128 * 1024, read::GzDecoder::new(file)))
    } else {
        Box::new(BufReader::with_capacity(128 * 1024, file))
    }
}

pub fn check_files_exist(args: &crate::Cli) {
    // check if the input files exist, and if not, exit gracefully

    if !args.bam.starts_with("s3") && !args.bam.starts_with("https://") {
        // TODO: We don't check for existence of remote files or their indexes
        // TODO: this could lead to nastier errors later, if the file happens to be unreachable
        // TODO: so maybe some time this should be properly implemented
        if !Path::new(&args.bam).exists() {
            error!("Alignment file not found: {}", args.bam);
            std::process::exit(1);
        }

        // Check for index file - support both .bai/.csi for BAM and .crai for CRAM
        if args.bam.ends_with(".cram") {
            let index = format!("{}.crai", args.bam);
            if !Path::new(&index).exists() {
                error!("Index file not found: {}", index);
                std::process::exit(1);
            }
        } else {
            // For BAM files, check for .bai first, then .csi as fallback
            let bai_index = format!("{}.bai", args.bam);
            let csi_index = format!("{}.csi", args.bam);
            if !Path::new(&bai_index).exists() && !Path::new(&csi_index).exists() {
                error!(
                    "Index file not found: {} or {}. Please create an index with 'samtools index'",
                    bai_index, csi_index
                );
                std::process::exit(1);
            }
        }
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
