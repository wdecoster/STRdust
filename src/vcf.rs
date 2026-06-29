use crate::Cli;
use crate::consensus::Consensus;
use crate::utils::{human_compare, levenshtein};
use log::debug;
use rust_htslib::faidx;
use std::cmp::Ordering;
use std::fmt;
use std::io::Read;

pub struct Allele {
    pub length: String, // length of the consensus sequence minus the length of the repeat sequence
    pub full_length: String, // length of the consensus sequence
    pub median_length: String, // median read length of the cluster, relative to the reference
    pub support: String, // number of reads supporting the allele
    pub std_dev: String, // standard deviation of the repeat length
    pub score: String,  // consensus score in the poa graph
    pub imprecise: bool, // cluster read-length spread exceeds the imprecision threshold
    pub seq: String,    // consensus sequence
}

impl Allele {
    pub fn from_consensus(consensus: Consensus, start: u32, end: u32) -> Allele {
        let ref_len = (end - start) as i32;
        match consensus.seq {
            Some(seq) => Allele {
                length: (seq.len() as i32 - ref_len).to_string(),
                full_length: seq.len().to_string(),
                median_length: (consensus.median_length as i32 - ref_len).to_string(),
                support: consensus.support.to_string(),
                std_dev: consensus.std_dev.to_string(),
                score: consensus.score.to_string(),
                imprecise: consensus.imprecise,
                seq,
            },
            None => Allele {
                length: ".".to_string(),
                full_length: ".".to_string(),
                median_length: ".".to_string(),
                support: consensus.support.to_string(),
                std_dev: ".".to_string(),
                score: ".".to_string(),
                imprecise: false,
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
    pub median_length: (String, String),
    pub support: (String, String),
    pub std_dev: (String, String),
    pub score: (String, String),
    pub somatic_info_field: String,
    pub outliers: String,
    pub n_clusters: String,
    pub dbscan_rb: String,
    pub time_taken: String,
    pub ps: Option<u32>, // phase set identifier
    pub flags: String,
    pub allele: (String, String),
    pub haploid: bool, // locus on a --haploid chromosome: report a single allele value (e.g. "1", not "1|1")
}

/// A chromosome listed (comma-separated) under `--haploid` is treated as haploid, so
/// genotypes at its loci are reported as a single allele value per the VCF specification.
pub fn chrom_is_haploid(args: &Cli, chrom: &str) -> bool {
    args.haploid
        .as_ref()
        .is_some_and(|h| h.split(',').any(|c| c == chrom))
}

impl VCFRecord {
    pub fn new(
        mut consenses: Vec<Consensus>,
        repeat_ref_sequence: String,
        all_insertions: Option<Vec<String>>,
        outlier_insertions: Option<Vec<String>>,
        n_clusters: Option<usize>,
        dbscan_rb: Option<String>,
        repeat: &crate::repeats::RepeatInterval,
        ps: Option<u32>,
        mut flag: Vec<String>,
        args: &Cli,
    ) -> VCFRecord {
        let haploid = chrom_is_haploid(args, &repeat.chrom);

        // On a haploid chromosome a single consensus is built, so the record carries one
        // allele; otherwise two consensus sequences are popped (in reverse, due to .pop()).
        let (allele1, allele2) = if haploid {
            let allele1 = Allele::from_consensus(
                consenses
                    .pop()
                    .expect("Failed getting allele1 from consenses"),
                repeat.start,
                repeat.end,
            );
            (allele1, None)
        } else {
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
            (allele1, Some(allele2))
        };

        let (genotype1, genotype2) = match &allele2 {
            Some(allele2) => {
                debug!(
                    "Genotyping {repeat}:{repeat_ref_sequence} with {} and {}",
                    allele1.seq, allele2.seq,
                );
                let gts = determine_genotypes(&repeat_ref_sequence, &allele1.seq, &allele2.seq);
                debug!(
                    "Genotype result for {repeat}: {}|{} (ref len: {}, allele1 len: {}, allele2 len: {})",
                    gts.0,
                    gts.1,
                    repeat_ref_sequence.len(),
                    allele1.seq.len(),
                    allele2.seq.len()
                );
                gts
            }
            None => {
                debug!("Genotyping haploid {repeat}:{repeat_ref_sequence} with {}", allele1.seq);
                let gt = determine_haploid_genotype(&repeat_ref_sequence, &allele1.seq);
                debug!("Haploid genotype result for {repeat}: {gt}");
                // The second slot is unused for haploid records (not emitted by Display).
                (gt, ".".to_string())
            }
        };

        let alts = match (genotype1.as_str(), genotype2.as_str()) {
            ("1", "0") | ("1", ".") | ("1", "1") => allele1.seq.clone(), // if both alleles are the same, only report one
            ("0", "1") | (".", "1") => allele2.as_ref().unwrap().seq.clone(),
            ("1", "2") => allele1.seq.clone() + "," + &allele2.as_ref().unwrap().seq,
            _ => ".".to_string(), // includes ./. and 0/0, and haploid 0 / .
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

        let n_clusters = match n_clusters {
            Some(n) => format!(";NCLUSTERS={n}"),
            None => "".to_string(),
        };

        let dbscan_rb = match dbscan_rb {
            Some(rb) => format!(";DBSCAN_RB={rb}"),
            None => "".to_string(),
        };

        // Flag the locus when either called allele's cluster has a wide read-length spread:
        // a single consensus length is then not representative (e.g. continuous/long-tailed loci).
        let allele2_imprecise = allele2.as_ref().map(|a| a.imprecise).unwrap_or(false);
        if (allele1.imprecise || allele2_imprecise) && !flag.iter().any(|f| f == "IMPRECISE_LENGTH")
        {
            flag.push("IMPRECISE_LENGTH".to_string());
        }
        let flags = if flag.is_empty() {
            "".to_string()
        } else {
            format!("{};", flag.join(";"))
        };
        let time_taken = repeat.time_field(args.debug);
        // The second slot of each FORMAT field holds the value of the second allele for diploid
        // records; for haploid records it stays "." and is not emitted by Display.
        let dot = || ".".to_string();
        VCFRecord {
            chrom: repeat.chrom.clone(),
            start: repeat.start,
            end: repeat.end,
            ref_seq: repeat_ref_sequence,
            alt_seq: Some(alts),
            length: (
                allele1.length,
                allele2
                    .as_ref()
                    .map(|a| a.length.clone())
                    .unwrap_or_else(&dot),
            ),
            full_length: (
                allele1.full_length,
                allele2
                    .as_ref()
                    .map(|a| a.full_length.clone())
                    .unwrap_or_else(&dot),
            ),
            median_length: (
                allele1.median_length,
                allele2
                    .as_ref()
                    .map(|a| a.median_length.clone())
                    .unwrap_or_else(&dot),
            ),
            support: (
                allele1.support,
                allele2
                    .as_ref()
                    .map(|a| a.support.clone())
                    .unwrap_or_else(&dot),
            ),
            std_dev: (
                allele1.std_dev,
                allele2
                    .as_ref()
                    .map(|a| a.std_dev.clone())
                    .unwrap_or_else(&dot),
            ),
            score: (
                allele1.score,
                allele2
                    .as_ref()
                    .map(|a| a.score.clone())
                    .unwrap_or_else(&dot),
            ),
            somatic_info_field,
            outliers,
            n_clusters,
            dbscan_rb,
            time_taken,
            ps,
            flags,
            allele: (genotype1.to_string(), genotype2.to_string()),
            haploid,
        }
    }

    pub fn quick_reference(
        repeat: &crate::repeats::RepeatInterval,
        repeat_ref_seq: &str,
        flags: Vec<String>,
        args: &crate::Cli,
    ) -> VCFRecord {
        let time_taken = repeat.time_field(args.debug);
        let flags_str = if flags.is_empty() {
            "".to_string()
        } else {
            format!("{};", flags.join(";"))
        };
        VCFRecord {
            chrom: repeat.chrom.clone(),
            start: repeat.start,
            end: repeat.end,
            ref_seq: repeat_ref_seq.to_string(),
            alt_seq: None,
            length: ("0".to_string(), "0".to_string()),
            full_length: (repeat_ref_seq.len().to_string(), repeat_ref_seq.len().to_string()),
            median_length: ("0".to_string(), "0".to_string()),
            support: (".".to_string(), ".".to_string()),
            std_dev: (".".to_string(), ".".to_string()),
            score: (".".to_string(), ".".to_string()),
            somatic_info_field: "".to_string(),
            outliers: "".to_string(),
            n_clusters: "".to_string(),
            dbscan_rb: "".to_string(),
            time_taken,
            ps: None,
            flags: flags_str,
            allele: ("0".to_string(), "0".to_string()),
            haploid: chrom_is_haploid(args, &repeat.chrom),
        }
    }

    pub fn missing_genotype(
        repeat: &crate::repeats::RepeatInterval,
        repeat_ref_seq: &str,
        support: String,
        args: &crate::Cli,
    ) -> VCFRecord {
        let time_taken = repeat.time_field(args.debug);
        VCFRecord {
            chrom: repeat.chrom.clone(),
            start: repeat.start,
            end: repeat.end,
            ref_seq: repeat_ref_seq.to_string(),
            alt_seq: Some(".".to_string()),
            length: (".".to_string(), ".".to_string()),
            full_length: (".".to_string(), ".".to_string()),
            median_length: (".".to_string(), ".".to_string()),
            support: (support, ".".to_string()),
            std_dev: (".".to_string(), ".".to_string()),
            score: (".".to_string(), ".".to_string()),
            somatic_info_field: "".to_string(),
            outliers: "".to_string(),
            n_clusters: "".to_string(),
            dbscan_rb: "".to_string(),
            time_taken,
            ps: None,
            flags: "".to_string(),
            allele: (".".to_string(), ".".to_string()),
            haploid: chrom_is_haploid(args, &repeat.chrom),
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
        let time_taken = repeat.time_field(args.debug);
        VCFRecord {
            chrom: repeat.chrom.clone(),
            start: repeat.start,
            end: repeat.end,
            ref_seq: repeat_ref_seq.to_string(),
            alt_seq: Some(seq.to_string()),
            length: ((seq.len() as u32 - (repeat.end - repeat.start)).to_string(), ".".to_string()),
            full_length: ((seq.len() as u32).to_string(), ".".to_string()),
            median_length: (
                (seq.len() as u32 - (repeat.end - repeat.start)).to_string(),
                ".".to_string(),
            ),
            support: ("1".to_string(), ".".to_string()),
            std_dev: (".".to_string(), ".".to_string()),
            score: (".".to_string(), ".".to_string()),
            somatic_info_field,
            outliers: "".to_string(),
            n_clusters: "".to_string(),
            dbscan_rb: "".to_string(),
            time_taken,
            ps,
            flags: if flag.is_empty() {
                "".to_string()
            } else {
                format!("{};", flag.join(";"))
            },
            allele: (".".to_string(), ".".to_string()),
            haploid: chrom_is_haploid(args, &repeat.chrom),
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

// Genotype a single allele on a haploid chromosome: "0" (reference), "1" (expansion/non-ref),
// or "." (missing). The VCF spec asks for a single allele value at haploid loci.
fn determine_haploid_genotype(repeat_ref_sequence: &str, allele: &str) -> String {
    if allele == "." {
        String::from(".")
    } else if is_similar_to_ref(allele, repeat_ref_sequence) {
        String::from("0")
    } else {
        String::from("1")
    }
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

impl VCFRecord {
    // GT field: a single allele value on haploid chromosomes ("1"), an allele pair otherwise ("1|0").
    fn genotype_field(&self) -> String {
        if self.haploid {
            self.allele.0.clone()
        } else {
            format!("{}|{}", self.allele.0, self.allele.1)
        }
    }

    // A per-allele FORMAT value: a single value on haploid chromosomes, a comma-separated pair otherwise.
    fn per_allele(&self, values: &(String, String)) -> String {
        if self.haploid {
            values.0.clone()
        } else {
            format!("{},{}", values.0, values.1)
        }
    }
}

impl fmt::Display for VCFRecord {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match &self.alt_seq {
            Some(alts) => {
                let (FORMAT, ps) = match self.ps {
                    Some(ps) => ("GT:RB:FRB:MRL:SUP:SC:PS", format!(":{}", ps)),
                    None => ("GT:RB:FRB:MRL:SUP:SC", "".to_string()),
                };
                write!(
                    f,
                    "{chrom}\t{start}\t.\t{ref}\t{alt}\t.\t.\t{flags}END={end};STDEV={stdev}{somatic}{outliers}{n_clusters}{dbscan_rb}{time_taken}\t{FORMAT}\t{gt}:{rb}:{frb}:{mrl}:{sup}:{sc}{ps}",
                    chrom = self.chrom,
                    start = self.start,
                    flags = self.flags,
                    end = self.end,
                    ref = self.ref_seq,
                    alt = alts,
                    stdev = self.per_allele(&self.std_dev),
                    somatic = self.somatic_info_field,
                    outliers = self.outliers,
                    n_clusters = self.n_clusters,
                    dbscan_rb = self.dbscan_rb,
                    time_taken = self.time_taken,
                    gt = self.genotype_field(),
                    rb = self.per_allele(&self.length),
                    frb = self.per_allele(&self.full_length),
                    mrl = self.per_allele(&self.median_length),
                    sup = self.per_allele(&self.support),
                    sc = self.per_allele(&self.score),
                )
            }
            None => {
                write!(
                    f,
                    "{chrom}\t{start}\t.\t{ref}\t.\t.\t.\t{flags}END={end};{somatic}\tGT:SUP\t{gt}:{sup}",
                    chrom = self.chrom,
                    start = self.start,
                    flags = self.flags,
                    end = self.end,
                    ref = self.ref_seq,
                    somatic = self.somatic_info_field,
                    gt = self.genotype_field(),
                    sup = self.per_allele(&self.support),
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
    println!(r#"##fileformat=VCFv4.5"#);
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
        r#"##INFO=<ID=STDEV,Number=.,Type=Integer,Description="Standard deviation of the repeat length per allele (one value on haploid chromosomes, two otherwise)">"#
    );
    println!(
        r#"##INFO=<ID=SEQS,Number=1,Type=String,Description="Sequences supporting the two alleles">"#
    );
    println!(
        r#"##INFO=<ID=OUTLIERS,Number=1,Type=String,Description="Outlier sequences much longer than the alleles">"#
    );
    println!(
        r#"##INFO=<ID=NCLUSTERS,Number=1,Type=Integer,Description="Number of read clusters found by the DBSCAN phasing strategy (>2 indicates a complex multi-population locus)">"#
    );
    println!(
        r#"##INFO=<ID=CLUSTERFAILURE,Number=0,Type=Flag,Description="If unphased input failed to cluster in two haplotype">"#
    );
    println!(
        r#"##INFO=<ID=EXPANSION_OUTLIER,Number=0,Type=Flag,Description="At least 2 reads are more than 2x longer than the longer called allele, suggesting a larger expansion missed by both alleles (only with --unphased)">"#
    );
    println!(
        r#"##INFO=<ID=IMPRECISE_LENGTH,Number=0,Type=Flag,Description="A called allele has a wide read-length spread (coefficient of variation > 0.2); the reported consensus length may not be representative">"#
    );
    println!(
        r#"##INFO=<ID=DISCORDANT_LENGTH,Number=0,Type=Flag,Description="With --phasing both: the reported Ward call and the DBSCAN call differ by more than 2x on the longer allele; see DBSCAN_RB (QC only, fires liberally)">"#
    );
    println!(
        r#"##INFO=<ID=DBSCAN_RB,Number=2,Type=Integer,Description="With --phasing both: DBSCAN allele lengths relative to reference, reported when DISCORDANT_LENGTH is set">"#
    );
    println!(
        r#"##INFO=<ID=QUICKREF,Number=0,Type=Flag,Description="Locus identified as homozygous reference via fast CIGAR check, alignment skipped">"#
    );
    if args.debug {
        println!(
            r#"##INFO=<ID=TIME,Number=1,Type=String,Description="Thread CPU time taken to genotype the repeat (debug mode only)">"#
        );
    }
    println!(r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#);
    println!(
        r#"##FORMAT=<ID=RB,Number=.,Type=Integer,Description="Repeat length per allele in bases relative to reference (one value on haploid chromosomes, two otherwise)">"#
    );
    println!(
        r#"##FORMAT=<ID=FRB,Number=.,Type=Integer,Description="Full repeat length per allele in bases (one value on haploid chromosomes, two otherwise)">"#
    );
    println!(
        r#"##FORMAT=<ID=MRL,Number=.,Type=Integer,Description="Median read length per allele's cluster in bases relative to reference, robust to a long length tail (one value on haploid chromosomes, two otherwise)">"#
    );
    println!(r#"##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set identifier">"#);
    println!(
        r#"##FORMAT=<ID=SUP,Number=.,Type=Integer,Description="Read support per allele (one value on haploid chromosomes, two otherwise)">"#
    );
    println!(
        r#"##FORMAT=<ID=SC,Number=.,Type=Integer,Description="Consensus score per allele (one value on haploid chromosomes, two otherwise)">"#
    );
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
        phasing_strategy: crate::PhasingStrategy::Ward,
        threads: 1,
        sample: None,
        haploid: Some(String::from("chr7")),
        debug: false,
        sorted: false,
        consensus_reads: 20,
        max_number_reads: 60,
        max_locus: None,
        alignment_all: false,
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
        phasing_strategy: crate::PhasingStrategy::Ward,
        threads: 1,
        sample: Some("test_sample".to_string()),
        haploid: Some(String::from("chr7")),
        debug: false,
        sorted: false,
        consensus_reads: 20,
        max_number_reads: 60,
        max_locus: None,
        alignment_all: false,
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

#[test]
fn test_determine_haploid_genotype() {
    let repeat_ref_sequence = "ATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATC";
    // Reference-matching allele -> 0
    assert_eq!(determine_haploid_genotype(repeat_ref_sequence, repeat_ref_sequence), "0");
    // Expanded allele -> 1
    let expanded = "ATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATC";
    assert_eq!(determine_haploid_genotype(repeat_ref_sequence, expanded), "1");
    // Missing allele -> .
    assert_eq!(determine_haploid_genotype(repeat_ref_sequence, "."), ".");
}

// A minimal record builder for testing the ploidy-aware Display output.
#[cfg(test)]
fn test_record(haploid: bool, allele: (&str, &str), alt: Option<&str>) -> VCFRecord {
    let pair = || ("10".to_string(), "20".to_string());
    VCFRecord {
        chrom: "chrX".to_string(),
        start: 100,
        end: 130,
        ref_seq: "ATCATCATCATCATCATCATCATCATCATC".to_string(),
        alt_seq: alt.map(|a| a.to_string()),
        length: pair(),
        full_length: pair(),
        median_length: pair(),
        support: pair(),
        std_dev: pair(),
        score: pair(),
        somatic_info_field: "".to_string(),
        outliers: "".to_string(),
        n_clusters: "".to_string(),
        dbscan_rb: "".to_string(),
        time_taken: "".to_string(),
        ps: None,
        flags: "".to_string(),
        allele: (allele.0.to_string(), allele.1.to_string()),
        haploid,
    }
}

#[test]
fn test_display_haploid_single_values() {
    // Haploid record: single GT value and single per-allele FORMAT values.
    let record = test_record(true, ("1", "."), Some("ATCATCATC"));
    let line = format!("{record}");
    let sample = line.split('\t').next_back().unwrap();
    assert_eq!(sample, "1:10:10:10:10:10");
    assert!(line.contains("STDEV=10\t") || line.contains("STDEV=10;"));
}

#[test]
fn test_display_diploid_pair_values() {
    // Diploid record: phased GT pair and comma-separated per-allele FORMAT values.
    let record = test_record(false, ("1", "0"), Some("ATCATCATC"));
    let line = format!("{record}");
    let sample = line.split('\t').next_back().unwrap();
    assert_eq!(sample, "1|0:10,20:10,20:10,20:10,20:10,20");
    assert!(line.contains("STDEV=10,20"));
}

#[test]
fn test_display_haploid_missing() {
    // Missing haploid call must be a single ".", not "./.".
    let record = test_record(true, (".", "."), Some("."));
    let line = format!("{record}");
    let sample = line.split('\t').next_back().unwrap();
    assert!(sample.starts_with("."), "expected single-dot GT, got {sample}");
    assert!(!sample.starts_with(".|."), "haploid missing must not be diploid ./.");
}
