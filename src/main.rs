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

/// What to optimise when deciding whether an allele counts as reference.
#[derive(clap::ValueEnum, Clone, Copy, Debug, PartialEq, Eq)]
pub enum Priority {
    /// Prioritise finding expanded alleles over resolving small ones. Loci whose reads all
    /// look near-reference are reported as reference without being genotyped, which is
    /// faster and costs resolution below ~3 bases. Long expansions are unaffected.
    Expanded,
    /// Report a variant unless the allele matches the reference exactly. Highest recall.
    Sensitive,
    /// Tolerate 2% of the reference length. Best balance of recall and precision.
    Balanced,
    /// Require a clear difference before reporting a variant. Highest precision.
    Precise,
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

    /// Minimal length of an insertion at the repeat junction for it to count towards the
    /// allele. Insertions shorter than this are read noise far more often than sequence:
    /// dropping to 1 in 02d9540 cost about 11 points of allele-length concordance against
    /// the GIAB HG002 truth set, and restoring a threshold recovers it
    #[arg(short, long, default_value_t = 3)]
    minlen: usize,

    /// minimal number of supporting reads per haplotype
    #[arg(short, long, default_value_t = 3)]
    support: usize,

    // ---------------------------------------------------------------------------------
    // Tuning knobs, hidden from --help on purpose. They exist so the benchmark can sweep a
    // parameter without a rebuild, which is what makes runs comparable by content rather
    // than by which branch someone had checked out. They are NOT a stable interface: once
    // the measurements settle, the winning value becomes the default and the flag goes.
    // Do not document these in the README and do not rely on them in scripts.
    // ---------------------------------------------------------------------------------
    // Sensitive-path only in practice: parse_cs is reached through find_insertions, so with
    // --mode fast this touches nothing except the clipped-read fallback - measured at 15
    // loci in 50,000 (0.03%), with every headline metric identical. In sensitive mode
    // narrowing it to 10 was worth about +0.5 on top of the other defaults, which is not
    // enough to justify a user-facing flag on a non-default path. Kept hidden because it
    // is the only handle on how much stray flank sequence parse_cs folds in, which is one
    // half of issue #24.
    //
    // NOTE an off-by-one that matters if this is ever narrowed further: the left flank of
    // the repeat-compressed reference is `flanking + 1` bases and the right is
    // `flanking - 1`, so the junction sits at offset flanking+1 while the window below is
    // centred on `flanking`. At +/-30 that is harmless; at 0 it excludes the junction
    // entirely and every locus no-calls.
    /// How far from the repeat/flank junction an insertion may sit and still be folded
    /// into the allele, in the sensitive path. Wider tolerates an aligner that places the
    /// insertion off the annotated boundary; narrower folds in less stray flank sequence
    #[arg(long, default_value_t = 30, hide = true)]
    junction_window: i32,

    /// What to prioritise when deciding whether an allele counts as reference.
    ///
    /// 'sensitive', 'balanced' and 'precise' change only the emitted genotype (GT): the
    /// measured allele lengths in RB, FRB and MRL are byte-identical under all three, so if
    /// you work from the lengths rather than the genotype they change nothing for you.
    ///
    /// 'expanded' is different in kind. It reports a locus as reference *without genotyping
    /// it* when every read looks near-reference, so at those loci RB and MRL are 0 rather
    /// than a measured value. Faster, and it gives up resolution below ~3 bases; long
    /// expansions are unaffected
    #[arg(long, value_enum, default_value_t = Priority::Balanced)]
    priority: Priority,

    // The three --ref-* knobs below are how --priority's three settings were derived, and
    // are kept hidden rather than deleted so the mapping can be re-derived if the consensus
    // changes again. Measured on 10,000 loci, fast + trim 0.35, recall/precision/F1 at
    // truth-variant loci:
    //     max_edits 0            99.0 / 54.2 / 70.1   -> priority sensitive
    //     divisor 50 (2% of REF) 85.5 / 58.7 / 69.6   -> priority balanced
    //     divisor 20, inclusive  39.9 / 91.8 / 55.7   -> priority precise
    //     divisor 20, strict     52.1 / 83.0 / 64.0   <- the old default, on no frontier
    // Note `inclusive` LOOSENS the rule (more alleles called reference), which is the
    // opposite of what the name suggests on first reading.
    /// Fixed edit distance an allele may differ from the reference and still be reported
    /// as reference. Negative defers to --priority
    #[arg(long, default_value_t = -1, hide = true, allow_hyphen_values = true)]
    ref_max_edits: i32,

    /// Divisor for the length-scaled reference-similarity threshold (20 = 5% of REF).
    /// 0 defers to --priority
    #[arg(long, default_value_t = 0, hide = true)]
    ref_edit_divisor: usize,

    /// Compare the edit distance with <= rather than <, so the stated tolerance is the
    /// real one rather than one edit tighter
    #[arg(long, default_value_t = false, hide = true)]
    ref_edit_inclusive: bool,

    /// Seed the POA graph with the read closest to the cluster median length instead of
    /// the first sampled read, whose indels would otherwise become the graph's backbone
    #[arg(long, default_value_t = true, hide = true)]
    poa_medoid_seed: bool,

    // 0.35 is the optimum of a six-point ladder (0.10/0.20/0.35/0.40/0.45/0.50) on 10,000
    // loci; it turns over on both sides, and above it the -1 bin grows, which is real
    // sequence being cut rather than one-read overhangs. The fraction maps to
    // ceil(fraction * n_reads), so its effective value is an integer that shifts with
    // depth - an absolute read count would be a cleaner parameterisation if this is ever
    // revisited.
    /// Fraction of a cluster's reads that must support an edge for the consensus to run
    /// through it at either end. Guards against rust-bio's consensus walking out along a
    /// single read's overhang. 0 keeps rust-bio's own endpoint choice
    #[arg(long, default_value_t = 0.35, hide = true)]
    poa_trim_fraction: f64,

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

    /// Reads are not phased: split them into haplotypes with the given strategy.
    /// 'ward': length-weighted Levenshtein with hierarchical clustering — the general
    /// choice, and the one to use when alleles differ mainly in length.
    /// 'dbscan': k-mer composition — separates alleles that differ in motif composition
    /// rather than in length, which Ward's length-weighted distance tends to merge.
    /// 'both': run Ward and additionally DBSCAN, reporting the Ward call but flagging
    /// substantial length discordance (DISCORDANT_LENGTH / DBSCAN_RB) — a QC mode for
    /// complicated regions where neither strategy should be trusted silently.
    ///
    /// A value is required: the strategies suit different situations and the right one is
    /// a property of the locus, not a default worth hiding.
    #[arg(long, value_name = "STRATEGY", value_enum)]
    unphased: Option<PhasingStrategy>,

    /// Identify poorly supported outlier expansions (only with --unphased)
    #[arg(long, default_value_t = false)]
    find_outliers: bool,

    /// Minimum fraction of reads required for a cluster to be considered a haplotype (only with --unphased)
    #[arg(long, default_value_t = 0.1)]
    min_haplotype_fraction: f32,

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
    /// 'fast': cut the repeat straight out of the alignment already in the BAM/CRAM
    /// (default).
    /// 'sensitive': re-align every read to a repeat-compressed reference.
    ///
    /// 'fast' is both quicker and more accurate against the GIAB HG002 truth set, including
    /// at long expansions, which was the case it was expected to lose: re-aligning
    /// reconstructs each read's allele by concatenating insertions near the junction, which
    /// folds in stray flanking sequence and cannot subtract deletions, while cutting from
    /// the alignment derives the allele from a reference span and gets both for free.
    #[arg(long, value_name = "MODE", value_enum, default_value_t = GenotypingMode::Fast)]
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

impl Cli {
    /// Net CIGAR length difference a read may show and still count as reference-like in the
    /// fast reference check (QUICKREF).
    ///
    /// Only `--priority expanded` loosens this. Measured on 50,000 loci: a tolerance of 3
    /// fires on 18,078 loci instead of 2,059 and cuts CPU by 15%, while *improving* exact
    /// concordance by 0.9 points and F1 by 2.5 - because a locus QUICKREF answers is one
    /// full genotyping does not get to answer wrongly. The cost is confined to the strata
    /// it should be: 1-10bp falls 57.4% -> 52.0% and 51-200bp 78.9% -> 73.3%, with long
    /// expansions untouched. Unlike the other three settings this changes *which loci are
    /// genotyped*, not just how the result is labelled.
    pub fn quickref_tolerance(&self) -> i64 {
        match self.priority {
            Priority::Expanded => 3,
            _ => 0,
        }
    }

    /// Whether reads must be split into haplotypes by clustering.
    pub fn is_unphased(&self) -> bool {
        self.unphased.is_some()
    }

    /// The clustering strategy, meaningful only when [`Cli::is_unphased`] is true.
    pub fn phasing_strategy(&self) -> PhasingStrategy {
        self.unphased.unwrap_or(PhasingStrategy::Ward)
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
    if args.find_outliers && !args.is_unphased() {
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
