use crate::Cli;
use crate::consensus::Consensus;
use distance::levenshtein;
use human_sort::compare as human_compare;
use log::debug;
use rust_htslib::faidx;
use std::cmp::Ordering;
use std::fmt;
use std::io::Read;

pub struct Allele {
    pub length: String, // length of the consensus sequence minus the length of the repeat sequence
    pub full_length: String, // length of the consensus sequence
    pub support: String, // number of reads supporting the allele
    pub std_dev: String, // standard deviation of the repeat length
    pub score: String,  // consensus score in the poa graph
    pub seq: String,    // consensus sequence
}

impl Allele {
    pub fn from_consensus(consensus: Consensus, start: u32, end: u32) -> Allele {
        match consensus.seq {
            Some(seq) => Allele {
                length: (seq.len() as i32 - ((end - start) as i32)).to_string(),
                full_length: seq.len().to_string(),
                support: consensus.support.to_string(),
                std_dev: consensus.std_dev.to_string(),
                score: consensus.score.to_string(),
                seq,
            },
            None => Allele {
                length: ".".to_string(),
                full_length: ".".to_string(),
                support: consensus.support.to_string(),
                std_dev: ".".to_string(),
                score: ".".to_string(),
                seq: ".".to_string(),
            },
        }
    }
}

pub struct VCFRecord {
    pub chrom: String,
    pub start: u32,
    pub end: u32,
    pub ref_seq: String,
    pub alt_seq: Option<String>,
    pub length: (String, String),
    pub full_length: (String, String),
    pub support: (String, String),
    pub std_dev: (String, String),
    pub score: (String, String),
    pub somatic_info_field: String,
    pub outliers: String,
    pub time_taken: String,
    pub ps: Option<u32>, // phase set identifier
    pub flags: String,
    pub allele: (String, String),
}

impl VCFRecord {
    pub fn new(
        mut consenses: Vec<Consensus>,
        repeat_ref_sequence: String,
        all_insertions: Option<Vec<String>>,
        outlier_insertions: Option<Vec<String>>,
        repeat: &crate::repeats::RepeatInterval,
        ps: Option<u32>,
        flag: Vec<String>,
        args: &Cli,
    ) -> VCFRecord {
        // since I use .pop() to format the two consensus sequences, the order is reversed
        let allele2 = Allele::from_consensus(
            consenses
                .pop()
                .expect("Failed getting allele2 from consenses"),
            repeat.start,
            repeat.end,
        );
        let allele1 = Allele::from_consensus(
            consenses
                .pop()
                .expect("Failed getting allele1 from consenses"),
            repeat.start,
            repeat.end,
        );

        debug!(
            "Genotyping {repeat}:{repeat_ref_sequence} with {} and {}",
            allele1.seq, allele2.seq,
        );
        let (genotype1, genotype2) =
            determine_genotypes(&repeat_ref_sequence, &allele1.seq, &allele2.seq);

        debug!(
            "Genotype result for {repeat}: {}|{} (ref len: {}, allele1 len: {}, allele2 len: {})",
            genotype1,
            genotype2,
            repeat_ref_sequence.len(),
            allele1.seq.len(),
            allele2.seq.len()
        );

        let alts = match (genotype1.as_str(), genotype2.as_str()) {
            ("1", "0") | ("1", ".") | ("1", "1") => allele1.seq, // if both alleles are the same, only report one
            ("0", "1") | (".", "1") => allele2.seq,
            ("1", "2") => allele1.seq + "," + &allele2.seq,
            _ => ".".to_string(), // includes ./. and 0/0
        };

        let somatic_info_field = match all_insertions {
            Some(somatic_insertions) => {
                format!(";SEQS={}", somatic_insertions.join(","))
            }
            None => "".to_string(),
        };

        let outliers = match outlier_insertions {
            Some(outlier_insertions) => {
                if outlier_insertions.is_empty() {
                    "".to_string()
                } else {
                    format!(";OUTLIERS={}", outlier_insertions.join(","))
                }
            }
            None => "".to_string(),
        };

        let flags = if flag.is_empty() {
            "".to_string()
        } else {
            format!("{};", flag.join(";"))
        };
        let time_taken = if args.debug {
            format!(
                ";TIME={}",
                chrono::Utc::now() - repeat.created.expect("Failed accessing timestamp")
            )
        } else {
            "".to_string()
        };
        VCFRecord {
            chrom: repeat.chrom.clone(),
            start: repeat.start,
            end: repeat.end,
            ref_seq: repeat_ref_sequence,
            alt_seq: Some(alts),
            length: (allele1.length, allele2.length),
            full_length: (allele1.full_length, allele2.full_length),
            support: (allele1.support, allele2.support),
            std_dev: (allele1.std_dev, allele2.std_dev),
            score: (allele1.score, allele2.score),
            somatic_info_field,
            outliers,
            time_taken,
            ps,
            flags,
            allele: (genotype1.to_string(), genotype2.to_string()),
        }
    }

    pub fn missing_genotype(
        repeat: &crate::repeats::RepeatInterval,
        repeat_ref_seq: &str,
        support: String,
        args: &crate::Cli,
    ) -> VCFRecord {
        let time_taken = if args.debug {
            format!(
                ";TIME={}",
                chrono::Utc::now() - repeat.created.expect("Failed accessing timestamp")
            )
        } else {
            "".to_string()
        };
        VCFRecord {
            chrom: repeat.chrom.clone(),
            start: repeat.start,
            end: repeat.end,
            ref_seq: repeat_ref_seq.to_string(),
            alt_seq: Some(".".to_string()),
            length: (".".to_string(), ".".to_string()),
            full_length: (".".to_string(), ".".to_string()),
            support: (support, ".".to_string()),
            std_dev: (".".to_string(), ".".to_string()),
            score: (".".to_string(), ".".to_string()),
            somatic_info_field: "".to_string(),
            outliers: "".to_string(),
            time_taken,
            ps: None,
            flags: "".to_string(),
            allele: (".".to_string(), ".".to_string()),
        }
    }

    pub fn single_read(
        seq: &str,
        repeat: &crate::repeats::RepeatInterval,
        repeat_ref_seq: &str,
        all_insertions: Option<Vec<String>>,
        ps: Option<u32>,
        flag: Vec<String>,
        args: &crate::Cli,
    ) -> VCFRecord {
        let somatic_info_field = match all_insertions {
            Some(somatic_insertions) => {
                format!(";SEQS={}", somatic_insertions.join(","))
            }
            None => "".to_string(),
        };
        let time_taken = if args.debug {
            format!(
                ";TIME={}",
                chrono::Utc::now() - repeat.created.expect("Failed accessing timestamp")
            )
        } else {
            "".to_string()
        };
        VCFRecord {
            chrom: repeat.chrom.clone(),
            start: repeat.start,
            end: repeat.end,
            ref_seq: repeat_ref_seq.to_string(),
            alt_seq: Some(seq.to_string()),
            length: ((seq.len() as u32 - (repeat.end - repeat.start)).to_string(), ".".to_string()),
            full_length: ((seq.len() as u32).to_string(), ".".to_string()),
            support: ("1".to_string(), ".".to_string()),
            std_dev: (".".to_string(), ".".to_string()),
            score: (".".to_string(), ".".to_string()),
            somatic_info_field,
            outliers: "".to_string(),
            time_taken,
            ps,
            flags: if flag.is_empty() {
                "".to_string()
            } else {
                format!("{};", flag.join(";"))
            },
            allele: (".".to_string(), ".".to_string()),
        }
    }
}

fn length_ratio(seq1: &str, seq2: &str) -> f32 {
    let length1 = seq1.len();
    let length2 = seq2.len();
    length1.min(length2) as f32 / length1.max(length2) as f32
}

fn determine_genotypes(
    repeat_ref_sequence: &str,
    allele1: &str,
    allele2: &str,
) -> (String, String) {
    // STR genotyping logic with clear decision matrix
    //
    // Strategy:
    // 1. Check if each allele is similar enough to reference (5% edit distance threshold)
    // 2. Check if alleles are similar to each other (for homozygous calls)
    // 3. Use match statement to handle all combinations explicitly
    //
    // Key principle: Both alleles can be called as reference (0|0) if they both
    // meet the similarity threshold, even if they differ slightly from each other

    // Handle missing alleles
    if allele1 == "." && allele2 == "." {
        return (String::from("."), String::from("."));
    } else if allele1 == "." {
        let gt2 = if is_similar_to_ref(allele2, repeat_ref_sequence) {
            "0"
        } else {
            "1"
        };
        return (String::from("."), String::from(gt2));
    } else if allele2 == "." {
        let gt1 = if is_similar_to_ref(allele1, repeat_ref_sequence) {
            "0"
        } else {
            "1"
        };
        return (String::from(gt1), String::from("."));
    }

    // Determine if each allele matches reference
    let allele1_is_ref = is_similar_to_ref(allele1, repeat_ref_sequence);
    let allele2_is_ref = is_similar_to_ref(allele2, repeat_ref_sequence);

    // Determine genotypes based on reference matching
    // Only compute alleles_same if needed (when neither matches ref)
    let genotypes = match (allele1_is_ref, allele2_is_ref) {
        // Both match reference -> 0|0
        (true, true) => ("0", "0"),

        // One matches ref, the other doesn't
        (true, false) => ("0", "1"),
        (false, true) => ("1", "0"),

        // Neither matches reference - now we need to check if they're similar to each other
        (false, false) => {
            if are_alleles_similar(allele1, allele2) {
                ("1", "1") // Similar to each other
            } else {
                ("1", "2") // Different from each other
            }
        }
    };

    (String::from(genotypes.0), String::from(genotypes.1))
}

// Check if an allele is similar enough to reference to be called as ref (0)
fn is_similar_to_ref(allele: &str, reference: &str) -> bool {
    // Exact match
    if allele == reference {
        debug!("  is_similar_to_ref: exact match");
        return true;
    }

    // Quick length-based filter: if length ratio <0.9, lengths too different
    // This is consistent with are_alleles_similar()
    let len_ratio = length_ratio(allele, reference);
    if len_ratio < 0.9 {
        debug!("  is_similar_to_ref: length ratio {} < 0.9, rejecting", len_ratio);
        return false;
    }

    // Use 5% edit distance threshold (same as original implementation)
    // Integer division means we're conservative
    // For 12bp: threshold = 0, so only exact matches pass
    // For 60bp: threshold = 3, so up to 2bp differences allowed
    // For 300bp: threshold = 15, so up to 14bp differences allowed
    let threshold = reference.len() / 20; // 5%, floor division

    // Optimization: if threshold is 0, we already checked for exact match above
    // No point computing expensive levenshtein distance - result will be false
    if threshold == 0 {
        debug!("  is_similar_to_ref: threshold is 0, rejecting");
        return false;
    }

    let edit_distance = levenshtein(allele, reference);
    let result = edit_distance < threshold; // Strictly less than (not <=)
    debug!(
        "  is_similar_to_ref: edit_distance {} vs threshold {}, result: {}",
        edit_distance, threshold, result
    );

    result
}

// Check if two alleles should be considered the same allele
// For STRs, focus on length similarity (same repeat count) rather than exact sequence
fn are_alleles_similar(allele1: &str, allele2: &str) -> bool {
    // Exact match
    if allele1 == allele2 {
        return true;
    }

    // For STRs, alleles with same length but different interruptions/SNPs
    // are often considered the same allele variant
    // Use length ratio >0.9 AND edit distance <10% as criteria
    let len_ratio = length_ratio(allele1, allele2);
    if len_ratio < 0.9 {
        return false;
    }

    // If lengths very similar, check edit distance
    // Allow up to 5% differences (accounts for sequencing errors, SNPs in repeat)
    let max_len = allele1.len().max(allele2.len());
    let threshold = std::cmp::max(1, max_len / 20); // 5%
    let edit_distance = levenshtein(allele1, allele2);

    edit_distance < threshold
}

impl fmt::Display for VCFRecord {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match &self.alt_seq {
            Some(alts) => {
                let (FORMAT, ps) = match self.ps {
                    Some(ps) => ("GT:RB:FRB:SUP:SC:PS", format!(":{}", ps)),
                    None => ("GT:RB:FRB:SUP:SC", "".to_string()),
                };
                write!(
                    f,
                    "{chrom}\t{start}\t.\t{ref}\t{alt}\t.\t.\t{flags}END={end};STDEV={sd1},{sd2}{somatic}{outliers}{time_taken}\t{FORMAT}\t{genotype1}|{genotype2}:{l1},{l2}:{fl1},{fl2}:{sup1},{sup2}:{score1},{score2}{ps}",
                    chrom = self.chrom,
                    start = self.start,
                    flags = self.flags,
                    end = self.end,
                    ref = self.ref_seq,
                    alt = alts,
                    l1 = self.length.0,
                    l2 = self.length.1,
                    fl1 = self.full_length.0,
                    fl2 = self.full_length.1,
                    sd1 = self.std_dev.0,
                    sd2 = self.std_dev.1,
                    somatic = self.somatic_info_field,
                    outliers = self.outliers,
                    time_taken = self.time_taken,
                    genotype1 = self.allele.0,
                    genotype2 = self.allele.1,
                    sup1 = self.support.0,
                    sup2 = self.support.1,
                    score1 = self.score.0,
                    score2 = self.score.1,
                )
            }
            None => {
                write!(
                    f,
                    "{chrom}\t{start}\t.\t{ref}\t.\t.\t.\tEND={end};{somatic}\tGT:SUP\t{genotype1}|{genotype2}:{sup1},{sup2}",
                    chrom = self.chrom,
                    start = self.start,
                    end = self.end,
                    ref = self.ref_seq,
                    somatic = self.somatic_info_field,
                    genotype1 = self.allele.0,
                    genotype2 = self.allele.1,
                    sup1 = self.support.0,
                    sup2 = self.support.1,
                )
            }
        }
    }
}

impl Ord for VCFRecord {
    fn cmp(&self, other: &Self) -> Ordering {
        human_compare(&self.chrom, &other.chrom).then(self.start.cmp(&other.start))
    }
}

impl PartialOrd for VCFRecord {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for VCFRecord {
    fn eq(&self, other: &Self) -> bool {
        (self.chrom.clone(), &self.start) == (other.chrom.clone(), &other.start)
    }
}

impl Eq for VCFRecord {}

pub fn write_vcf_header(args: &Cli) {
    println!(r#"##fileformat=VCFv4.2"#);
    // get absolute path to fasta file
    let path = std::fs::canonicalize(&args.fasta)
        .unwrap_or_else(|err| panic!("Failed getting absolute path to fasta: {err}"));
    println!(r#"##reference={}"#, path.to_str().expect("Failed converting path to string"));
    // get the version of this crate
    let version = env!("CARGO_PKG_VERSION");
    println!(r#"##source=STRdust v{}"#, version);
    // write all arguments to the header
    println!(r#"##command=STRdust {}"#, std::env::args().collect::<Vec<String>>().join(" "));
    // call faidx to make sure the fasta index exists, we'll need this anyway when genotyping
    let _ = faidx::Reader::from_path(&args.fasta)
        .unwrap_or_else(|err| panic!("Failed opening fasta: {err}"));

    let mut fai_file =
        std::fs::File::open(format!("{0}.fai", args.fasta)).expect("Can't open .fai file");
    // parse the fasta index file
    let mut buf = String::new();
    fai_file
        .read_to_string(&mut buf)
        .expect("Can't read fai file");
    for contig in buf.lines() {
        let mut contig = contig.split_whitespace();
        let name = contig.next().unwrap();
        let length = contig.next().unwrap().parse::<usize>().unwrap();
        println!(r#"##contig=<ID={},length={}>"#, name, length);
    }
    println!(
        r#"##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the repeat interval">"#
    );
    println!(
        r#"##INFO=<ID=STDEV,Number=2,Type=Integer,Description="Standard deviation of the repeat length">"#
    );
    println!(
        r#"##INFO=<ID=SEQS,Number=1,Type=String,Description="Sequences supporting the two alleles">"#
    );
    println!(
        r#"##INFO=<ID=OUTLIERS,Number=1,Type=String,Description="Outlier sequences much longer than the alleles">"#
    );
    println!(
        r#"##INFO=<ID=CLUSTERFAILURE,Number=0,Type=Flag,Description="If unphased input failed to cluster in two haplotype">"#
    );
    if args.debug {
        println!(
            r#"##INFO=<ID=TIME,Number=1,Type=String,Description="Time taken to genotype the repeat">"#
        );
    }
    println!(r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#);
    println!(
        r#"##FORMAT=<ID=RB,Number=2,Type=Integer,Description="Repeat length of the two alleles in bases relative to reference">"#
    );
    println!(
        r#"##FORMAT=<ID=FRB,Number=2,Type=Integer,Description="Full repeat length of the two alleles in bases">"#
    );
    println!(r#"##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set identifier">"#);
    println!(r#"##FORMAT=<ID=SUP,Number=2,Type=Integer,Description="Read support per allele">"#);
    println!(r#"##FORMAT=<ID=SC,Number=2,Type=Integer,Description="Consensus score per allele">"#);
    let name = match &args.sample {
        Some(name) => name,
        None => {
            // use basename of bam and remove file extension

            std::path::Path::new(&args.bam)
                .file_stem()
                .unwrap()
                .to_str()
                .unwrap()
        }
    };
    println!("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{name}",);
}

#[cfg(test)]
#[test]
fn test_write_vcf_header_from_bam() {
    let args = Cli {
        bam: String::from("test_data/small-test-phased.bam"),
        fasta: String::from("test_data/chr7.fa.gz"),
        region: None,
        region_file: None,
        pathogenic: false,
        minlen: 5,
        support: 1,
        somatic: false,
        unphased: true,
        find_outliers: false,
        min_haplotype_fraction: 0.1,
        threads: 1,
        sample: None,
        haploid: Some(String::from("chr7")),
        debug: false,
        sorted: false,
        consensus_reads: 20,
        max_number_reads: 60,
        max_locus: None,
    };
    write_vcf_header(&args);
}

#[test]
fn test_write_vcf_header_from_name() {
    let args = Cli {
        bam: String::from("test_data/small-test-phased.bam"),
        fasta: String::from("test_data/chr7.fa.gz"),
        region: None,
        region_file: None,
        pathogenic: false,
        minlen: 5,
        support: 1,
        somatic: false,
        unphased: true,
        find_outliers: false,
        min_haplotype_fraction: 0.1,
        threads: 1,
        sample: Some("test_sample".to_string()),
        haploid: Some(String::from("chr7")),
        debug: false,
        sorted: false,
        consensus_reads: 20,
        max_number_reads: 60,
        max_locus: None,
    };
    write_vcf_header(&args);
}

#[test]
fn test_determine_genotypes() {
    let repeat_ref_sequence = "ATCATCATCATC";
    let allele1 = "ATCATCATCATC";
    let allele2 = "ATCATCATCATC";
    let (genotype1, genotype2) = determine_genotypes(repeat_ref_sequence, allele1, allele2);
    assert_eq!(genotype1, "0");
    assert_eq!(genotype2, "0");
}

#[test]
fn test_determine_genotypes2() {
    let repeat_ref_sequence = "ATCATCATCATC";
    let allele1 = "ATCATCATCATC";
    let allele2 = "ATCATCATCATG";
    let (genotype1, genotype2) = determine_genotypes(repeat_ref_sequence, allele1, allele2);
    assert_eq!(genotype1, "0");
    assert_eq!(genotype2, "1");
}

#[test]
fn test_determine_genotypes3() {
    let repeat_ref_sequence = "ATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATC";
    let allele1 = "ATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATC";
    let allele2 = "ATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATG";
    let (genotype1, genotype2) = determine_genotypes(repeat_ref_sequence, allele1, allele2);
    assert_eq!(genotype1, "1");
    assert_eq!(genotype2, "0");
}

#[test]
fn test_determine_genotypes4() {
    let repeat_ref_sequence = "GGAGGAGGAGGAGGA";
    let allele1 = "GGAGGAGGAGGAGGAGGAGGAGGA";
    let allele2 = "GGAGGAGGAGGAGGAGGAGGAGGAAGGAGGAGGAGGAGGAGGAGGA";
    let (genotype1, genotype2) = determine_genotypes(repeat_ref_sequence, allele1, allele2);
    assert_eq!(genotype1, "1");
    assert_eq!(genotype2, "2");
}

#[test]
fn test_determine_genotypes5() {
    let repeat_ref_sequence = "GGAGGAGGAGGAGGA";
    let allele1 = "GGAGGAGGAGGAGGAGGAGGAGGAAGGAGGAGGAGGAGGAGGAGGA";
    let allele2 = "GGAGGAGGAGGAGGAGGAGGAGGAAGGAGGAGGAGGAGGAGGAGGA";
    let (genotype1, genotype2) = determine_genotypes(repeat_ref_sequence, allele1, allele2);
    assert_eq!(genotype1, "1");
    assert_eq!(genotype2, "1");
}

#[test]
fn test_determine_genotypes6() {
    let repeat_ref_sequence = "GGAGGAGGAGGAGGA";
    let allele1 = "GGAGGAGGAGGAGGAGGAGGAGGAAGGAGGAGGAGGAGGAGGAGGA";
    let allele2 = "GGAGGAGGAGGAGGAGGAGGAGGAAGGAGGAGGAGGAGGAGGAGTA";
    let (genotype1, genotype2) = determine_genotypes(repeat_ref_sequence, allele1, allele2);
    assert_eq!(genotype1, "1");
    assert_eq!(genotype2, "1");
}

#[test]
fn test_determine_genotypes7() {
    let repeat_ref_sequence = "GGAGGAGGAGGAGGA";
    let allele1 = "GGAGGAGGAGGAGGAGGAGGAGGAAGGAGGAGGAGGAGGAGGAGGA";
    let allele2 = ".";
    let (genotype1, genotype2) = determine_genotypes(repeat_ref_sequence, allele1, allele2);
    assert_eq!(genotype1, "1");
    assert_eq!(genotype2, ".");
}
