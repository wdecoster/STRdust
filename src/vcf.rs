use crate::consensus::Consensus;
use distance::levenshtein;
use human_sort::compare as human_compare;
use log::debug;
use rust_htslib::faidx;
use std::cmp::Ordering;
use std::fmt;
use std::io::Read;
use crate::Cli;


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
            length: (
                (seq.len() as u32 - (repeat.end - repeat.start)).to_string(),
                ".".to_string(),
            ),
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
    // if the consensus is very similar to the reference the variant is considered ref
    // for this I use a threshold of 5% of the length of the repeat sequence in the reference
    // e.g. if the repeat is 300bp in the reference this will allow an edit distance of 15
    // not sure if these numbers require further tuning
    // note that is an integer division, i.e. floor division
    // the next_alt variable dictates which genotype code the genotype2 can be in the case it is not the same as genotype1
    // this avoids a 0/2 genotype
    // the levenshtein distance should not be calculated if there is a large difference in length between the sequences, as this is computationally expensive and unnecessary
    // the ratio is used to determine if the sequences are similar enough to be considered the same allele
    let (genotype1, next_alt) = if allele1 == "." {
        (".", "1")
    // if the length ratio is less than 0.9, there is a large difference in length between the sequences
    // therefore the allele is considered different from the reference
    } else if length_ratio(allele1, repeat_ref_sequence) < 0.9 {
        ("1", "2")
    // if the length ratio is between 0.9 and 1 (with the latter indicating the same length) then we check the edit distance
    // if the edit distance is less than 5% of the length of the repeat sequence in the reference, the allele is considered the same as the reference
    } else if repeat_ref_sequence == allele1
        || levenshtein(allele1, repeat_ref_sequence) < repeat_ref_sequence.len() / 20
    {
        ("0", "1")
    // if the edit distance is larger than the threshold, the allele is considered different from the reference
    } else {
        ("1", "2")
    };

    let genotype2 = if allele2 == "." {
        "."
    // first check if the alt alleles are the same. if so, they have the same genotype
    } else if allele2 == allele1 {
        genotype1
    // if the length ratio is less than 0.9, there is a large difference between the sequence and they cannot be the same allele
    // but then we still have to check if the allele is the same as the reference
    } else if length_ratio(allele1, allele2) < 0.9 {
        // if the allele is different from the reference and from the first allele, we pick the next_alt value
        if length_ratio(allele2, repeat_ref_sequence) < 0.9 {
            next_alt
        // if the allele is similar to the reference, we pick the 0 value
        } else if levenshtein(allele2, repeat_ref_sequence) < repeat_ref_sequence.len() / 20 {
            "0"
        // if the allele is different from the reference, we pick the next_alt value
        } else {
            next_alt
        }
    // if the length_ratio of allele1 and allele2 was not less than 0.9, we check the edit distance
    // if the edit distance is less that 5%, we considered the alleles the same
    } else if levenshtein(allele2, allele1) < allele1.len() / 20 {
        genotype1
    // if the edit distance is larger than the threshold, the allele is compared with the reference
    } else if length_ratio(allele2, repeat_ref_sequence) < 0.9 {
        next_alt
    // if the allele is not too different in length from the reference, we check its edit distance
    // the threshold is the same as for the first allele, 5% of the length of the repeat sequence in the reference
    } else if levenshtein(allele2, repeat_ref_sequence) < repeat_ref_sequence.len() / 20 {
        "0"
    // if the allele is actualy different from the reference, we pick the next_alt value
    } else {
        next_alt
    };
    (String::from(genotype1), String::from(genotype2))
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
    println!(
        r#"##reference={}"#,
        path.to_str().expect("Failed converting path to string")
    );
    // get the version of this crate
    let version = env!("CARGO_PKG_VERSION");
    println!(r#"##source=STRdust v{}"#, version);
    // write all arguments to the header
    println!(
        r#"##command=STRdust {}"#,
        std::env::args().collect::<Vec<String>>().join(" ")
    );
    // call faidx to make sure the fasta index exists, we'll need this anyway when genotyping
    let _ =
        faidx::Reader::from_path(&args.fasta).unwrap_or_else(|err| panic!("Failed opening fasta: {err}"));
    
    let mut fai_file = std::fs::File::open(format!("{0}.fai", args.fasta)).expect("Can't open .fai file");
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
        println!(r#"##INFO=<ID=TIME,Number=1,Type=String,Description="Time taken to genotype the repeat">"#);
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
            let name = std::path::Path::new(&args.bam)
                .file_stem()
                .unwrap()
                .to_str()
                .unwrap();
            name
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
            threads: 1,
            sample: None,
            haploid: Some(String::from("chr7")),
            debug: false,
            sorted: false,
            consensus_reads: 20,
            max_number_reads: 60,
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
            threads: 1,
            sample: Some("test_sample".to_string()),
            haploid: Some(String::from("chr7")),
            debug: false,
            sorted: false,
            consensus_reads: 20,
            max_number_reads: 60,
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
