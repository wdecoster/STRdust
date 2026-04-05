use flate2::read;
use log::error;
use std::cmp;
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

/// Levenshtein edit distance between two strings.
/// Operates on bytes directly since STR sequences are ASCII.
pub fn levenshtein(a: &str, b: &str) -> usize {
    let a = a.as_bytes();
    let b = b.as_bytes();
    let b_len = b.len();
    let mut prev: Vec<usize> = (0..=b_len).collect();
    let mut curr = vec![0; b_len + 1];

    for (i, &ca) in a.iter().enumerate() {
        curr[0] = i + 1;
        for (j, &cb) in b.iter().enumerate() {
            let cost = usize::from(ca != cb);
            curr[j + 1] = cmp::min(cmp::min(prev[j + 1] + 1, curr[j] + 1), prev[j] + cost);
        }
        std::mem::swap(&mut prev, &mut curr);
    }
    prev[b_len]
}

/// Compare two strings in natural (human) order, where embedded numbers
/// are compared numerically rather than lexicographically.
pub fn human_compare(a: &str, b: &str) -> cmp::Ordering {
    let mut a_chars = a.chars().peekable();
    let mut b_chars = b.chars().peekable();

    loop {
        match (a_chars.peek(), b_chars.peek()) {
            (None, None) => return cmp::Ordering::Equal,
            (None, Some(_)) => return cmp::Ordering::Less,
            (Some(_), None) => return cmp::Ordering::Greater,
            (Some(&ca), Some(&cb)) => {
                if ca.is_ascii_digit() && cb.is_ascii_digit() {
                    // Extract full numeric segments
                    let na = extract_number(&mut a_chars);
                    let nb = extract_number(&mut b_chars);
                    match na.cmp(&nb) {
                        cmp::Ordering::Equal => continue,
                        ord => return ord,
                    }
                }
                match ca.cmp(&cb) {
                    cmp::Ordering::Equal => {
                        a_chars.next();
                        b_chars.next();
                    }
                    ord => return ord,
                }
            }
        }
    }
}

fn extract_number(chars: &mut std::iter::Peekable<std::str::Chars>) -> u64 {
    let mut n: u64 = 0;
    while let Some(&c) = chars.peek() {
        if c.is_ascii_digit() {
            n = n.saturating_mul(10).saturating_add(c as u64 - '0' as u64);
            chars.next();
        } else {
            break;
        }
    }
    n
}

/// Cross-platform cache directory.
pub fn cache_dir() -> std::path::PathBuf {
    std::env::var_os("XDG_CACHE_HOME")
        .map(std::path::PathBuf::from)
        .or_else(|| {
            #[cfg(target_os = "macos")]
            {
                std::env::var_os("HOME").map(|h| std::path::PathBuf::from(h).join("Library/Caches"))
            }
            #[cfg(target_os = "windows")]
            {
                std::env::var_os("LOCALAPPDATA").map(std::path::PathBuf::from)
            }
            #[cfg(not(any(target_os = "macos", target_os = "windows")))]
            {
                std::env::var_os("HOME").map(|h| std::path::PathBuf::from(h).join(".cache"))
            }
        })
        .unwrap_or_else(std::env::temp_dir)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_levenshtein_identical() {
        assert_eq!(levenshtein("ACGT", "ACGT"), 0);
    }

    #[test]
    fn test_levenshtein_empty_strings() {
        assert_eq!(levenshtein("", ""), 0);
    }

    #[test]
    fn test_levenshtein_one_empty() {
        assert_eq!(levenshtein("ACGT", ""), 4);
        assert_eq!(levenshtein("", "ACGT"), 4);
    }

    #[test]
    fn test_levenshtein_single_substitution() {
        assert_eq!(levenshtein("ACGT", "ACTT"), 1);
    }

    #[test]
    fn test_levenshtein_single_insertion() {
        assert_eq!(levenshtein("ACGT", "ACGGT"), 1);
    }

    #[test]
    fn test_levenshtein_single_deletion() {
        assert_eq!(levenshtein("ACGT", "ACT"), 1);
    }

    #[test]
    fn test_levenshtein_completely_different() {
        assert_eq!(levenshtein("AAAA", "TTTT"), 4);
    }

    #[test]
    fn test_levenshtein_symmetric() {
        assert_eq!(levenshtein("ACGT", "TGCA"), levenshtein("TGCA", "ACGT"));
    }

    #[test]
    fn test_levenshtein_repeated_motif() {
        assert_eq!(levenshtein("CAGCAGCAG", "CAGCAGCAGCAG"), 3);
    }

    #[test]
    fn test_levenshtein_classic_example() {
        assert_eq!(levenshtein("kitten", "sitting"), 3);
    }

    #[test]
    fn test_levenshtein_single_char() {
        assert_eq!(levenshtein("A", "T"), 1);
        assert_eq!(levenshtein("A", "A"), 0);
    }
}
