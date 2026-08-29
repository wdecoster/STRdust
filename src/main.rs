#![allow(non_snake_case)]
use clap::{Parser, ValueEnum};
use log::{info, warn};
use std::path::PathBuf;

pub mod bam_pool;
pub mod batching;
pub mod call;
pub mod consensus;
pub mod dbscan;
pub mod features;
pub mod genotype;
pub mod motif;
pub mod parse_bam;
pub mod phase_insertions;
pub mod repeats;
pub mod utils;
pub mod vcf;

/// How the repeat sequence of a read is recovered.
#[derive(ValueEnum, Clone, Copy, Debug, PartialEq, Eq)]
pub enum GenotypingMode {
    /// Re-align every read to a reference with the repeat excised and take the sequence
    /// that fails to align. Slow, but recovers alleles from reads that the original
    /// alignment clipped or placed badly, which is what large expansions look like.
    Sensitive,
    /// Cut the repeat straight out of the CIGAR of the alignment in the BAM/CRAM. Orders
    /// of magnitude faster, but only reads that span the locus can contribute; loci left
    /// without enough spanning reads fall back to `sensitive`.
    Fast,
}

/// Strategy for splitting unphased reads into haplotypes.
#[derive(ValueEnum, Clone, Copy, Debug, PartialEq, Eq)]
pub enum PhasingStrategy {
    /// Length-weighted Levenshtein distance + Ward hierarchical clustering.
    Ward,
    /// k-mer composition feature vectors + DBSCAN.
    Dbscan,
    /// QC mode: report the Ward call but additionally run DBSCAN and flag substantial
    /// length discordance (DISCORDANT_LENGTH / DBSCAN_RB) for review.
    Both,
}

// The arguments end up in the Cli struct
#[derive(Parser, Debug)]
#[command(author, version, about = "Tool to genotype STRs from long reads", long_about = None)]
pub struct Cli {
    /// reference genome
    #[arg(value_parser = is_file)]
    fasta: String,

    /// BAM or CRAM file to call STRs in
    #[arg(value_parser = is_file)]
    bam: String,

    /// Region string to genotype expansion in (format: chr:start-end)
    #[arg(short, long)]
    region: Option<String>,

    /// Bed file with region(s) to genotype expansion(s) in  
    #[arg(short = 'R', long, value_parser = is_file)]
    region_file: Option<String>,

    /// Genotype the pathogenic STRs from STRchive
    #[arg(long, default_value_t = false)]
    pathogenic: bool,

    /// minimal length of insertion/deletion operation
    #[arg(short, long, default_value_t = 1)]
    minlen: usize,

    /// minimal number of supporting reads per haplotype
    #[arg(short, long, default_value_t = 3)]
    support: usize,

    /// Minimum mapping quality of a read to be used. Lower it (down to 0) to keep
    /// ambiguously mapped reads, which matters in segmental duplications
    #[arg(long, default_value_t = 10)]
    mapq: u8,

    /// Number of parallel threads to use
    #[arg(short, long, default_value_t = 1)]
    threads: usize,

    /// Sample name to use in VCF header, if not provided, the bam file name is used
    #[arg(long)]
    sample: Option<String>,

    /// Print information on somatic variability
    #[arg(long, default_value_t = false)]
    somatic: bool,

    /// Reads are not phased
    #[arg(long, default_value_t = false)]
    unphased: bool,

    /// Identify poorly supported outlier expansions (only with --unphased)
    #[arg(long, default_value_t = false)]
    find_outliers: bool,

    /// Minimum fraction of reads required for a cluster to be considered a haplotype (only with --unphased)
    #[arg(long, default_value_t = 0.1)]
    min_haplotype_fraction: f32,

    /// Strategy for splitting unphased reads into haplotypes (only with --unphased).
    /// 'ward': length-weighted Levenshtein + hierarchical clustering (default).
    /// 'dbscan': k-mer composition features + DBSCAN (experimental, robust to length-variable expansions).
    /// 'both': QC mode that reports the Ward call but also runs DBSCAN and flags substantial
    /// length discordance (DISCORDANT_LENGTH / DBSCAN_RB) for review.
    #[arg(long = "phasing", value_name = "STRATEGY", value_enum, default_value_t = PhasingStrategy::Ward)]
    phasing_strategy: PhasingStrategy,

    /// comma-separated list of haploid (sex) chromosomes
    #[arg(long)]
    haploid: Option<String>,

    /// Debug mode
    #[arg(long, default_value_t = false)]
    debug: bool,

    /// Sort output by chrom, start and end
    #[arg(long, default_value_t = false)]
    sorted: bool,

    /// Max number of reads to use to generate consensus alt sequence
    #[arg(long, default_value_t = 20)]
    consensus_reads: usize,

    /// Max number of reads to extract per locus from the bam file for genotyping (use -1 for all reads)
    #[arg(long, default_value_t = 60, allow_hyphen_values = true)]
    max_number_reads: isize,

    /// Maximum locus size to consider (intervals larger than this will be filtered out)
    #[arg(long)]
    max_locus: Option<u32>,

    /// Always use full alignment (disable fast reference check via CIGAR)
    #[arg(long, default_value_t = false)]
    alignment_all: bool,

    /// How to recover the repeat sequence from a read.
    /// 'sensitive': re-align every read to a repeat-compressed reference (default).
    /// 'fast': cut the repeat straight out of the alignment already in the BAM/CRAM.
    #[arg(long, value_name = "MODE", value_enum, default_value_t = GenotypingMode::Sensitive)]
    mode: GenotypingMode,

    /// How far outside the annotated interval an insertion may sit and still count towards
    /// the allele with --mode fast. Aligners place a large insertion inconsistently, often
    /// tens of bases off the repeat; the default mirrors the tolerance of the sensitive path
    #[arg(long, default_value_t = 20)]
    fast_flank: u32,
}

impl Cli {
    /// Whether the repeat sequence is cut out of the existing alignment rather than
    /// recovered by re-aligning the read.
    pub fn is_fast_mode(&self) -> bool {
        self.mode == GenotypingMode::Fast
    }
}

fn is_file(pathname: &str) -> Result<String, String> {
    if pathname.starts_with("http")
        || pathname.starts_with("https://")
        || pathname.starts_with("s3")
    {
        return Ok(pathname.to_string());
    }

    let path = PathBuf::from(pathname);
    if path.is_file() {
        Ok(pathname.to_string())
    } else {
        Err(format!("Input file {} is invalid", path.display()))
    }
}

fn main() {
    env_logger::init();
    let args = Cli::parse();
    if args.find_outliers && !args.unphased {
        warn!("--find-outliers is only effective with --unphased");
    }
    if args.haploid.is_some() {
        warn!(
            "As of v0.20.0, genotypes on --haploid chromosomes are reported as a single allele \
             value (e.g. GT '1', or '.' when missing) per the VCF specification, instead of the \
             previous diploid representation ('1/1', './.'). Per-allele FORMAT/INFO fields (RB, \
             FRB, MRL, SUP, SC, STDEV) likewise carry a single value at these loci. Downstream \
             tools that assumed diploid genotypes may need updating."
        );
    }
    if args.is_fast_mode() && args.minlen != 1 {
        warn!(
            "--minlen has no effect with --mode fast: the allele is read off the alignment \
             rather than collected from insertion operations, so there is no indel length to \
             filter on."
        );
    }
    info!("Collected arguments: {args:?}");
    call::genotype_repeats(args);
}

#[cfg(test)]
#[ctor::ctor(unsafe)]
fn init() {
    env_logger::init();
}

#[test]
fn verify_app() {
    use clap::CommandFactory;
    Cli::command().debug_assert()
}
