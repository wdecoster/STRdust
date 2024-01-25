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
    ) -> VCFRecord {
        // since I use .pop() to format the two consensus sequences, the order is reversed
        let allele2 = Allele::from_consensus(consenses.pop().unwrap(), repeat.start, repeat.end);
        let allele1 = Allele::from_consensus(consenses.pop().unwrap(), repeat.start, repeat.end);

        debug!(
            "Genotyping {repeat}:{repeat_ref_sequence} with {} and {}",
            allele1.seq, allele2.seq,
        );
        // if the consensus is very similar to the reference the variant is considered ref
        // for this I use a threshold of 5% of the length of the repeat sequence in the reference
        // e.g. if the repeat is 300bp in the reference this will allow an edit distance of 15
        // not sure if these numbers require further tuning
        // note that is an integer division, i.e. floor division
        // the next_alt variable dictates which genotype code the genotype2 can be in the case it is not the same as genotype1
        // this avoids a 0/2 genotype
        let (genotype1, next_alt) = if allele1.seq == "." {
            (".", "1")
        } else if levenshtein(&allele1.seq, &repeat_ref_sequence) < repeat_ref_sequence.len() / 20 {
            ("0", "1")
        } else {
            ("1", "2")
        };

        let genotype2 = if allele2.seq == "." {
            "."
        } else if levenshtein(&allele2.seq, &repeat_ref_sequence) < repeat_ref_sequence.len() / 20 {
            "0"
        } else if levenshtein(&allele2.seq, &allele1.seq) < allele1.seq.len() / 20 {
            genotype1
        } else {
            next_alt
        };

        let alts = match (genotype1, genotype2) {
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
            ps,
            flags,
            allele: (genotype1.to_string(), genotype2.to_string()),
        }
    }

    pub fn missing_genotype(
        repeat: &crate::repeats::RepeatInterval,
        repeat_ref_seq: &str,
        support: String,
    ) -> VCFRecord {
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
            ps: None,
            flags: "".to_string(),
            allele: (".".to_string(), ".".to_string()),
        }
    }
}

impl fmt::Display for VCFRecord {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match &self.alt_seq {
            Some(alts) => {
                let (FORMAT, ps) = match self.ps {
                    Some(ps) => ("GT:SUP:SC:PS", format!(":{}", ps)),
                    None => ("GT:SUP:SC", "".to_string()),
                };
                write!(
                    f,
                    "{chrom}\t{start}\t.\t{ref}\t{alt}\t.\t.\t{flags}END={end};RB={l1},{l2};FRB={fl1},{fl2};STDEV={sd1},{sd2}{somatic}{outliers}\t{FORMAT}\t{genotype1}|{genotype2}:{sup1},{sup2}:{score1},{score2}{ps}",
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

pub fn write_vcf_header(fasta: &str, bam: &str, sample: &Option<String>) {
    println!(r#"##fileformat=VCFv4.2"#);
    // get absolute path to fasta file
    let path = std::fs::canonicalize(fasta)
        .unwrap_or_else(|err| panic!("Failed getting absolute path to fasta: {err}"));
    println!(
        r#"##reference={}"#,
        path.to_str().expect("Failed converting path to string")
    );
    // get the version of this crate
    let version = env!("CARGO_PKG_VERSION");
    println!(r#"##source=STRdust v{}"#, version);
    // call faidx to make sure the fasta index exists, we'll need this anyway when genotyping
    let _ =
        faidx::Reader::from_path(fasta).unwrap_or_else(|err| panic!("Failed opening fasta: {err}"));

    let mut fai_file = std::fs::File::open(format!("{fasta}.fai")).expect("Can't open file");
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
        r#"##INFO=<ID=RB,Number=2,Type=Integer,Description="Repeat length of the two alleles in bases relative to reference">"#
    );
    println!(
        r#"##INFO=<ID=FRB,Number=2,Type=Integer,Description="Full repeat length of the two alleles in bases">"#
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
    println!(r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#);
    println!(r#"##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set identifier">"#);
    println!(r#"##FORMAT=<ID=SUP,Number=2,Type=Integer,Description="Read support per allele">"#);
    println!(r#"##FORMAT=<ID=SC,Number=2,Type=Integer,Description="Consensus score per allele">"#);
    let name = match sample {
        Some(name) => name,
        None => {
            // use basename of bam and remove file extension
            let name = std::path::Path::new(&bam)
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
    write_vcf_header(
        "test_data/chr7.fa.gz",
        "test_data/small-test-phased.bam",
        &None,
    );
}

#[test]
fn test_write_vcf_header_from_name() {
    write_vcf_header(
        "test_data/chr7.fa.gz",
        "test_data/small-test-phased.bam",
        &Some("test_sample".to_string()),
    );
}
