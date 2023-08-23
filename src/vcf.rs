use crate::consensus::Consensus;
use distance::levenshtein;
use human_sort::compare as human_compare;
use log::debug;
use rust_htslib::faidx;
use std::cmp::Ordering;
use std::fmt;
use std::io::Read;

pub struct VCFRecord {
    pub chrom: String,
    pub start: u32,
    pub end: u32,
    pub ref_seq: String,
    pub alt_seq: Option<String>,
    pub length: (String, String),
    pub support: (String, String),
    pub std_dev: (String, String),
    pub score: (String, String),
    pub somatic_info_field: String,
    pub flags: String,
    pub allele: (String, String),
}

impl VCFRecord {
    pub fn new(
        mut consenses: Vec<Consensus>,
        repeat_ref_sequence: String,
        somatic: bool,
        all_insertions: Vec<String>,
        repeat: crate::repeats::RepeatInterval,
        flag: Vec<String>,
    ) -> VCFRecord {
        // since I use .pop() to format the two consensus sequences, the order is reversed
        let (length2, alt2, support2, std_dev2, score2) = consenses
            .pop()
            .unwrap()
            .format_lengths(repeat.start, repeat.end);
        let (length1, alt1, support1, std_dev1, score1) = consenses
            .pop()
            .unwrap()
            .format_lengths(repeat.start, repeat.end);
        debug!("Genotyping {repeat}:{repeat_ref_sequence} with {alt1} and {alt2}");
        let allele1 = if alt1 == "." {
            "."
        // if the consensus is very similar to the reference the variant is considered ref
        // for this I use a threshold of 5% of the length of the repeat sequence in the reference
        // e.g. if the repeat is 300bp in the reference this will allow an edit distance of 15
        // not sure if these numbers require further tuning
        // note that is an integer division, i.e. floor division
        } else if levenshtein(&alt1, &repeat_ref_sequence) < repeat_ref_sequence.len() / 20 {
            "0"
        } else {
            "1"
        };

        let allele2 = if alt2 == "." {
            "."
        } else if levenshtein(&alt2, &repeat_ref_sequence) < repeat_ref_sequence.len() / 20 {
            "0"
        } else if alt2 == alt1 {
            "1"
        } else {
            "2"
        };

        let alts = match (allele1, allele2) {
            ("1", "0") | ("1", ".") | ("1", "1") => alt1, // if both alleles are the same, only report one
            ("0", "1") | (".", "1") => alt2,
            ("1", "2") => alt1 + "," + &alt2,
            _ => ".".to_string(), // includes ./. and 0/0
        };

        let somatic_info_field = if somatic {
            format!(";SEQS={}", all_insertions.join("|"))
        } else {
            "".to_string()
        };

        let flags = flag.join(";");
        VCFRecord {
            chrom: repeat.chrom,
            start: repeat.start,
            end: repeat.end,
            ref_seq: repeat_ref_sequence,
            alt_seq: Some(alts),
            length: (length1, length2),
            support: (support1, support2),
            std_dev: (std_dev1, std_dev2),
            score: (score1, score2),
            somatic_info_field,
            flags,
            allele: (allele1.to_string(), allele2.to_string()),
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
            support: (support, ".".to_string()),
            std_dev: (".".to_string(), ".".to_string()),
            score: (".".to_string(), ".".to_string()),
            somatic_info_field: "".to_string(),
            flags: "".to_string(),
            allele: (".".to_string(), ".".to_string()),
        }
    }
}

impl fmt::Display for VCFRecord {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match &self.alt_seq {
            Some(alts) => {
                write!(
            f,
            "{chrom}\t{start}\t.\t{ref_seq}\t{alt_seq}\t.\t.\t\
            {flags}END={end};RL={length1}|{length2};SUPP={support1}|{support2};STDEV={std_dev1}|{std_dev2};CONSENSUS_SCORE={score1}|{score2};{somatic_info_field}\
            \tGT\t{allele1}|{allele2}",
            chrom = self.chrom,
            start = self.start,
            flags = self.flags,
            end = self.end,
            ref_seq = self.ref_seq,
            alt_seq = alts,
            length1 = self.length.0,
            length2 = self.length.1,
            support1 = self.support.0,
            support2 = self.support.1,
            std_dev1 = self.std_dev.0,
            std_dev2 = self.std_dev.1,
            score1 = self.score.0,
            score2 = self.score.1,
            somatic_info_field = self.somatic_info_field,
            allele1 = self.allele.0,
            allele2 = self.allele.1,
        )
            }
            None => {
                write!(
                    f,
                    "{chrom}\t{start}\t.\t{ref_seq}\t.\t.\t.\tEND={end};SUPP={support1}|{support2};{somatic_info_field}\tGT\t{allele1}|{allele2}",
                    chrom = self.chrom,
                    start = self.start,
                    end = self.end,
                    ref_seq = self.ref_seq,
                    support1 = self.support.0,
                    support2 = self.support.1,
                    somatic_info_field = self.somatic_info_field,
                    allele1 = self.allele.0,
                    allele2 = self.allele.1,
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

pub fn write_vcf_header(fasta: &str, bam: &str, sample: Option<String>) {
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
        r#"##INFO=<ID=RL,Number=2,Type=Integer,Description="Repeat length of the two alleles">"#
    );
    println!(
        r#"##INFO=<ID=SUPP,Number=2,Type=Integer,Description="Number of reads supporting the two alleles">"#
    );
    println!(
        r#"##INFO=<ID=STDEV,Number=2,Type=Integer,Description="Standard deviation of the repeat length">"#
    );
    println!(
        r#"##INFO=<ID=SEQS,Number=1,Type=String,Description="Sequences supporting the two alleles">"#
    );
    println!(r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#);
    println!(r#"##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set identifier">"#);
    let name = match sample {
        Some(name) => name,
        None => {
            // use basename of bam and remove file extension
            let name = std::path::Path::new(&bam)
                .file_stem()
                .unwrap()
                .to_str()
                .unwrap();
            name.to_string()
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
        None,
    );
}

#[test]
fn test_write_vcf_header_from_name() {
    write_vcf_header(
        "test_data/chr7.fa.gz",
        "test_data/small-test-phased.bam",
        Some("test_sample".to_string()),
    );
}
