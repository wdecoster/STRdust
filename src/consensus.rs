use log::{debug, warn};
use poasta::align::PoastaAligner;
use poasta::align::astar::AstarEngine;
use poasta::align::cost_models::linear::Linear;
use poasta::align::traits::{AlignResult, AlignmentEngine};
use poasta::graph::consensus::{heaviest_bundle_consensus, trim_consensus_ends};
use poasta::graph::poa::{POAGraph, POANodeIndex};
use poasta::graph::traits::{GraphBase, GraphWithNodeOrdering};
use rand::SeedableRng;
use rand::rngs::StdRng;
use rand::seq::IteratorRandom;
use std::fmt;

/// Fixed seed for the read downsampling RNG. Downsampling is purely a
/// performance/memory measure, so the choice of subset should not make the
/// genotype non-reproducible: seeding with a constant makes every run (and
/// every comparison between runs) deterministic and independent of locus or
/// thread ordering.
const DOWNSAMPLE_SEED: u64 = 42;

/// A called allele is flagged IMPRECISE_LENGTH when the coefficient of variation
/// (std_dev / mean) of its cluster's read lengths exceeds this threshold. 0.2 (20%
/// spread) separates clean, near-fixed-length alleles from continuous/long-tailed
/// length distributions (e.g. the GOLGA8A repeat) where a single consensus length
/// is not a faithful representation of the underlying reads. Hardcoded for now;
/// can be promoted to a CLI argument later if cohort experience warrants it.
const IMPRECISE_LENGTH_CV: f64 = 0.2;

/// `PoaTies::MedianLength` only arbitrates alleles whose `score` consensus lies within this
/// many bases of the reference length. Reads through long expansions run short of the true
/// allele, so there the median read length pulled calls too short: on the long-expansion
/// truth set picking by median cost a point of exact matches against `score`, while
/// gating at 20 kept its gain on 10,000 random loci (20, 50 and 100 performed alike).
const MEDIAN_LENGTH_MAX_DEVIATION: usize = 20;

/// Scoring for the POA (poasta) that builds an allele consensus.
///
/// poasta works with costs rather than scores: a match is free and mismatches and gap
/// bases are penalised, with a linear gap model. `Default` is the plain algorithm: the old
/// rust-bio scoring (match 3, mismatch -4, linear gap -12) converted to costs and doubled to
/// stay integral (mismatch = 2 * (3 + 4) = 14, gap = 2 * (3 / 2 + 12) = 27), no medoid seed,
/// no trimming and poasta's own tie-break. The command line defaults differ: cheaper costs
/// (3/2), which simulated repeat reads favoured at high error rates and for long alleles and
/// the truth set confirmed, and [`PoaTies::MedianLength`]. Affine gaps did clearly worse on
/// repeats in simulation.
#[derive(Clone, Copy, Debug)]
pub struct PoaScoring {
    pub mismatch: u8,
    pub gap: u8,
    /// Seed the POA graph with the read whose length is closest to the cluster median,
    /// rather than whichever read happened to be sampled first.
    pub medoid_seed: bool,
    /// Fraction of the cluster's reads that must cover the first and last base of the
    /// consensus. 0 disables trimming.
    pub trim_fraction: f64,
    /// How the heaviest bundle resolves two equally heavy in-edges.
    pub ties: PoaTies,
    /// Reads outside mean +/- this many standard deviations are dropped before the POA.
    pub outlier_std_factor: f64,
    /// Length outliers are only dropped when the standard deviation is at least this.
    pub outlier_min_std: usize,
}

/// How the heaviest bundle resolves two equally heavy in-edges.
///
/// A tie is common in homopolymers, where reads of different lengths split over insertion
/// nodes. On 10,000 HG002 truth-set loci, `score` made alleles too long more often than too
/// short (1144 vs 675 in fast mode), `coverage` the reverse (695 vs 1124).
#[derive(clap::ValueEnum, Clone, Copy, Debug, PartialEq, Eq)]
pub enum PoaTies {
    /// poasta's own order: the predecessor with the higher path score, i.e. the longer path
    Score,
    /// The predecessor covered by more reads, see [`heaviest_bundle_coverage_ties`]
    Coverage,
    /// Build both and keep the consensus whose length is closest to the cluster's median
    /// read length, `score` on a tie. Only near the reference length, see
    /// [`MEDIAN_LENGTH_MAX_DEVIATION`]
    MedianLength,
}

impl Default for PoaScoring {
    fn default() -> Self {
        Self {
            mismatch: 14,
            gap: 27,
            medoid_seed: false,
            trim_fraction: 0.0,
            ties: PoaTies::Score,
            outlier_std_factor: 2.0,
            outlier_min_std: 5,
        }
    }
}

/// Index of the read whose length is closest to the cluster median.
///
/// The first read added to a POA graph is its backbone: every other read is aligned onto it,
/// so its indels are structurally privileged and end up in the consensus. Taking whichever
/// read was sampled first makes that an arbitrary choice. The read at the median length is
/// the one least likely to drag the consensus off the cluster's centre. Ties go to the
/// earlier read, so the choice stays deterministic.
fn medoid_index(seqs: &[Vec<u8>]) -> usize {
    let mut lengths: Vec<usize> = seqs.iter().map(|s| s.len()).collect();
    lengths.sort_unstable();
    let median = lengths[lengths.len() / 2];
    seqs.iter()
        .enumerate()
        .min_by_key(|(_, s)| s.len().abs_diff(median))
        .map(|(i, _)| i)
        .unwrap_or(0)
}

/// Heaviest bundle consensus as poasta's `heaviest_bundle_consensus`, except for how ties
/// on edge weight are broken.
///
/// poasta resolves a tie by the higher accumulated path score, which is the longer path, so
/// in a homopolymer the consensus systematically runs long. E.g. 4 reads with 16 A's and
/// 2 + 2 reads with 17 and 18 A's that merge into the same last node give the G after the
/// run two in-edges of weight 8, and the 18 A path wins. Here a tie goes to the predecessor
/// covered by more reads (the 16th A, covered by all 8), then falls back to poasta's order.
fn heaviest_bundle_coverage_ties(graph: &POAGraph<u32>) -> Vec<POANodeIndex<u32>> {
    if graph.is_empty() {
        return Vec::new();
    }
    let (start, end) = (graph.start_node(), graph.end_node());
    // edges from the start and to the end sentinel weigh 0, so take the larger side
    let coverage = |v: POANodeIndex<u32>| {
        let into: usize = graph.in_edges(v).map(|(_, e)| e.weight).sum();
        let out: usize = graph.out_edges(v).map(|(_, e)| e.weight).sum();
        into.max(out)
    };
    let n = graph.node_count();
    let mut score = vec![0usize; n];
    let mut pred: Vec<Option<POANodeIndex<u32>>> = vec![None; n];
    for rank in 0..n {
        let v = graph.rank_to_node(rank);
        if v == start {
            continue;
        }
        // (edge weight, predecessor coverage, predecessor score, lowest predecessor rank)
        let mut best = None;
        for (u, edge) in graph.in_edges(v) {
            let u_rank = graph.node_rank(u);
            let key = (edge.weight, coverage(u), score[u_rank], std::cmp::Reverse(u_rank));
            if best.is_none_or(|b| key > b) {
                best = Some(key);
            }
        }
        if let Some((w, _, s, std::cmp::Reverse(u_rank))) = best {
            score[rank] = s + w;
            pred[rank] = Some(graph.rank_to_node(u_rank));
        }
    }
    let mut path = Vec::new();
    let mut curr = end;
    while let Some(p) = pred[graph.node_rank(curr)] {
        path.push(curr);
        curr = p;
    }
    path.reverse();
    path.retain(|v| *v != start && *v != end);
    path
}

/// Build the POA graph from `seqs` (the read at `seed` first) and return its consensus
/// together with the cost of aligning that consensus back to the graph.
///
/// The consensus is the heaviest bundle, trimmed at both ends to bases covered by at least
/// `trim_fraction` of the reads. Without the trimming the path always runs to the deepest
/// sink, so a single read extending one base past the others makes the allele one base too
/// long (reported and fixed upstream in poasta#37). If trimming leaves nothing, the
/// untrimmed consensus is used rather than an empty allele.
fn poa_consensus(
    seqs: &[Vec<u8>],
    seed: usize,
    median_length: usize,
    ref_length: usize,
    poa: PoaScoring,
) -> Result<(Vec<u8>, i32), String> {
    let costs = Linear::new(poa.mismatch, poa.gap);
    let mut aligner =
        PoastaAligner::new(AstarEngine::<POAGraph<u32>, _>::new(costs), POAGraph::<u32>::new());
    let order = std::iter::once(&seqs[seed])
        .chain(
            seqs.iter()
                .enumerate()
                .filter(|(i, _)| *i != seed)
                .map(|(_, s)| s),
        )
        .map(|s| s.as_slice());
    aligner.align_all(order).map_err(|e| e.to_string())?;
    let graph = aligner.graph();

    let spell = |path: Vec<POANodeIndex<u32>>| -> Vec<u8> {
        let trimmed = trim_consensus_ends(graph, &path, poa.trim_fraction);
        let path = if trimmed.is_empty() {
            &path[..]
        } else {
            trimmed
        };
        path.iter().map(|&n| graph.node_symbol(n)).collect()
    };
    let consensus = match poa.ties {
        PoaTies::Score => spell(heaviest_bundle_consensus(graph)),
        PoaTies::Coverage => spell(heaviest_bundle_coverage_ties(graph)),
        PoaTies::MedianLength => {
            let by_score = spell(heaviest_bundle_consensus(graph));
            if by_score.len().abs_diff(ref_length) > MEDIAN_LENGTH_MAX_DEVIATION {
                by_score
            } else {
                let by_coverage = spell(heaviest_bundle_coverage_ties(graph));
                if by_coverage.len().abs_diff(median_length)
                    < by_score.len().abs_diff(median_length)
                {
                    by_coverage
                } else {
                    by_score
                }
            }
        }
    };

    let cost = AstarEngine::<POAGraph<u32>, _>::new(costs)
        .align(graph, &consensus)
        .map_err(|e| e.to_string())?
        .alignment_score();
    Ok((consensus, cost as i32))
}

#[derive(Clone)]
pub struct Consensus {
    pub seq: Option<String>,
    pub support: usize,
    pub std_dev: usize,
    /// Cost of aligning the consensus back to its POA graph (lower is better); 0 for a
    /// single-read consensus and -1 when there is no consensus.
    pub score: i32,
    /// Median full length of the cluster reads (before length-outlier removal),
    /// reported as MRL. More robust to a long tail than the POA consensus length.
    pub median_length: usize,
    /// True when the cluster's read-length coefficient of variation exceeds
    /// [`IMPRECISE_LENGTH_CV`]; surfaced as the IMPRECISE_LENGTH locus flag.
    pub imprecise: bool,
}

impl Default for Consensus {
    fn default() -> Consensus {
        Consensus {
            seq: None,
            support: 0,
            std_dev: 0,
            score: -1,
            median_length: 0,
            imprecise: false,
        }
    }
}

impl fmt::Display for Consensus {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match &self.seq {
            Some(seq) => write!(
                f,
                "seq: {}, support: {}, std_dev: {}, score: {}",
                seq, self.support, self.std_dev, self.score
            ),
            None => write!(f, "seq: None, support: 0, std_dev: 0, score: -1"),
        }
    }
}

pub fn consensus(
    seqs: &[String],
    support: usize,
    consensus_reads: usize,
    repeat: &crate::repeats::RepeatInterval,
    poa: PoaScoring,
) -> Consensus {
    if seqs.is_empty() {
        return Consensus::default();
    }
    let num_reads_ = seqs.len();
    // Length statistics over all cluster reads (before length-outlier removal), so the
    // reported median and the imprecision verdict reflect the true read distribution.
    let (mean, std_dev, median_length) = length_stats(seqs);
    let imprecise = mean > 0 && (std_dev as f64 / mean as f64) > IMPRECISE_LENGTH_CV;
    let seqs = remove_outliers(seqs, mean, std_dev, &poa, repeat);
    let num_reads = seqs.len();
    debug!("{repeat}: Kept {}/{} reads after dropping outliers", num_reads, num_reads_);
    if num_reads < support {
        Consensus { seq: None, support: num_reads, std_dev, score: -1, median_length, imprecise }
    } else if consensus_reads == 1 {
        // if only on read should be used to generate the consensus, the consensus is a randomly selected read
        let seq = seqs
            .into_iter()
            .choose(&mut StdRng::seed_from_u64(DOWNSAMPLE_SEED))
            .unwrap();
        Consensus {
            seq: Some(seq.to_string()),
            support: num_reads,
            std_dev,
            score: 0,
            median_length,
            imprecise,
        }
    } else {
        // if there are more than <consensus_reads> reads, downsample before taking the consensus
        // for performance and memory reasons
        let seqs_bytes = if num_reads > consensus_reads {
            debug!("{repeat}: Too many reads, downsampling to {consensus_reads}");
            seqs.into_iter()
                .sample(&mut StdRng::seed_from_u64(DOWNSAMPLE_SEED), consensus_reads)
                .into_iter()
                .map(|seq| seq.bytes().collect::<Vec<u8>>())
                .collect::<Vec<Vec<u8>>>()
        } else {
            seqs.iter()
                .map(|seq| seq.bytes().collect::<Vec<u8>>())
                .collect::<Vec<Vec<u8>>>()
        };
        log::info!("Creating consensus for {repeat}");
        let seed = if poa.medoid_seed {
            medoid_index(&seqs_bytes)
        } else {
            0
        };
        debug!("{repeat}: seeding POA graph with read {seed} of {}", seqs_bytes.len());
        // the annotated repeat spans end - start + 1 reference bases
        let ref_length = (repeat.end - repeat.start + 1) as usize;
        let (consensus, score) =
            match poa_consensus(&seqs_bytes, seed, median_length, ref_length, poa) {
                Ok(result) => result,
                Err(e) => {
                    // never lose the allele over an aligner failure: the seed read is the
                    // best single-read stand-in, and -1 marks the score as unavailable
                    warn!("{repeat}: POA failed ({e}), using the seed read as consensus");
                    (seqs_bytes[seed].clone(), -1)
                }
            };
        debug!("Created consensus");
        Consensus {
            seq: Some(std::str::from_utf8(&consensus).unwrap().to_string()),
            support: num_reads,
            std_dev,
            score,
            median_length,
            imprecise,
        }
    }
}

/// Mean, (floored) standard deviation, and median of the sequence lengths.
fn length_stats(seqs: &[String]) -> (usize, usize, usize) {
    let mut lengths = seqs.iter().map(|x| x.len()).collect::<Vec<usize>>();
    let mean = lengths.iter().sum::<usize>() / lengths.len();
    let variance = lengths
        .iter()
        .map(|x| (*x as isize - mean as isize).pow(2) as usize)
        .sum::<usize>()
        / lengths.len();
    // casting to usize floors the std_dev to the integer below, rather than rounding to nearest,
    // so this keeps the std_dev slightly smaller than it really is, but not by a lot
    let std_dev = (variance as f64).sqrt() as usize;
    lengths.sort_unstable();
    let median = if lengths.len() % 2 == 0 {
        (lengths[lengths.len() / 2] + lengths[lengths.len() / 2 - 1]) / 2
    } else {
        lengths[lengths.len() / 2]
    };
    (mean, std_dev, median)
}

fn remove_outliers<'a>(
    seqs: &'a [String],
    mean: usize,
    std_dev: usize,
    poa: &PoaScoring,
    repeat: &crate::repeats::RepeatInterval,
) -> Vec<&'a String> {
    // remove sequences that are shorter or longer than <factor> standard deviations from the
    // mean, except if the stdev is small
    debug!("{repeat}: mean: {}, std_dev: {}", mean, std_dev);
    if std_dev < poa.outlier_min_std {
        debug!("std_dev < {}, not removing any outliers", poa.outlier_min_std);
        seqs.iter().collect::<Vec<&String>>()
    } else {
        let window = (poa.outlier_std_factor * std_dev as f64) as usize;
        // avoid underflowing usize
        let min_val = mean.saturating_sub(window);
        let max_val = mean + window;
        debug!("Removing outliers outside of [{},{}]", min_val, max_val);
        seqs.iter()
            .filter(|seq| seq.len() > min_val && seq.len() < max_val)
            .collect::<Vec<&String>>()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dummy_repeat() -> crate::repeats::RepeatInterval {
        crate::repeats::RepeatInterval {
            chrom: "chr1".to_string(),
            start: 1,
            end: 100,
            created: None,
        }
    }

    #[test]
    fn test_consensus_reports_median_and_not_imprecise_when_tight() {
        // near-uniform lengths (~30 bp): low CV -> not imprecise, median ~30
        let seqs: Vec<String> = (0..10).map(|i| "A".repeat(30 + i % 2)).collect();
        let cons = consensus(&seqs, 2, 1, &dummy_repeat(), PoaScoring::default());
        assert!(!cons.imprecise, "tight length distribution should not be imprecise");
        assert!((29..=31).contains(&cons.median_length), "median {}", cons.median_length);
    }

    #[test]
    fn test_consensus_flags_imprecise_when_length_spread_is_wide() {
        // continuous/long-tailed lengths: CV well above 0.2 -> imprecise
        let seqs: Vec<String> = vec![30, 50, 120, 300, 900, 2700]
            .into_iter()
            .map(|l| "A".repeat(l))
            .collect();
        let cons = consensus(&seqs, 2, 1, &dummy_repeat(), PoaScoring::default());
        assert!(cons.imprecise, "wide length distribution should be flagged imprecise");
    }

    #[test]
    fn test_consensus_ignores_single_read_overhangs() {
        // 1 of 15 reads extends one base past the others, at either end: the untrimmed
        // heaviest bundle runs out along it, the trimmed consensus does not
        let base = "CAG".repeat(20);
        let trim = PoaScoring { trim_fraction: 0.35, ..PoaScoring::default() };
        for overhang in [format!("{base}C"), format!("T{base}")] {
            let mut seqs = vec![base.clone(); 14];
            seqs.push(overhang);
            let cons = consensus(&seqs, 1, 20, &dummy_repeat(), trim);
            assert_eq!(cons.seq.as_deref(), Some(base.as_str()));
            let untrimmed = consensus(&seqs, 1, 20, &dummy_repeat(), PoaScoring::default());
            assert_eq!(untrimmed.seq.map(|s| s.len()), Some(61));
        }
    }

    #[test]
    fn test_consensus_homopolymer_ties() {
        // HG002 reads at chr1:74729637 (truth: reference, 16 A's): 4 reads with 16 A's,
        // 2 with 17 and 2 with 18. Score tie-breaking calls 18 A's, coverage gives 16, and
        // so does the median length (50 bases) choosing between the two, as long as the
        // allele is near the reference length
        let (a16, a17, a18) = ("A".repeat(16), "A".repeat(17), "A".repeat(18));
        let read = |a: &str| format!("CCCTGTCTC{a}GCAGCAGACTAGAGAAGTCATTAAA");
        let seqs: Vec<String> = [&a17, &a16, &a18, &a16, &a18, &a17, &a16, &a16]
            .iter()
            .map(|a| read(a))
            .collect();
        let poa = PoaScoring { medoid_seed: true, trim_fraction: 0.35, ..PoaScoring::default() };
        let by_score = consensus(&seqs, 1, 20, &dummy_repeat(), poa);
        assert_eq!(by_score.seq, Some(read(&a18)));
        let by_coverage =
            consensus(&seqs, 1, 20, &dummy_repeat(), PoaScoring { ties: PoaTies::Coverage, ..poa });
        assert_eq!(by_coverage.seq, Some(read(&a16)));
        let median = PoaScoring { ties: PoaTies::MedianLength, ..poa };
        let reference = crate::repeats::RepeatInterval { end: 50, ..dummy_repeat() };
        let by_median = consensus(&seqs, 1, 20, &reference, median);
        assert_eq!(by_median.seq, Some(read(&a16)));
        // 50 bases off the 100 base reference, the score consensus is kept
        let far_from_reference = consensus(&seqs, 1, 20, &dummy_repeat(), median);
        assert_eq!(far_from_reference.seq, Some(read(&a18)));
    }

    #[test]
    fn test_consensus() {
        // I created this test because these sequences segfaulted on bianca
        let seqs = vec![
            "CAGACAGGCAGGCAGGCAGGCAGGCAGGCAGGCAGACAGGCAGCCAGGCAGGCAGGCAGG".to_string(),
            "CAGGCAGGCAGGCAGGCAGGCAGGCAGGCAGACAGGCAGCCAGGCAGGCAGGCAGG".to_string(),
            "AGACAGACAGGCAGGCAGGCAGGCAGGCAGGCAGGCAGACAGGCAGCCAGGCAGGCAGGCAGGC".to_string(),
            "GGCAGGCAGGCAGGCAGGCAGGCAGACAGGCAGCCAGGCAGGCAGGCAGGCAGACAGAAG".to_string(),
            "CAGGCAGGCAGGCAGGCAGGCAGGCAGGCAGACAGGCAGCCAGGCAGGCAGGCAGG".to_string(),
            "CAGGCAGGCAGGCAGGCAGGCAGGCAGGCAGACAGGCAGCCAGGCAGGCAGGCAGG".to_string(),
            "CAGACAGGCAGGCAGGCAGGCAGGCAGGCAGGCAGACAGGCAGCCAGGCAGGCAGGCAGG".to_string(),
            "CAGGCAGGCAGGCAGGCAGGCAGGCAGGCAGACAGGCAGCCAGGCAGGCAGGCAGG".to_string(),
            "AGACAGACAGGCAGGCAGGCAGGCAGGCAGGCAGGCAGACAGGCAGCCAGGCAGGCAGGCAGGC".to_string(),
            "CAGGCAGGCAGGCAGGCAGGCAGGCAGGCAGACAGGCAGCCAGGCAGGCAGGCAGG".to_string(),
            "CAGGCAGGCAGGCAGGCAGGCAGGCAGGCAGACAGGCAGCCAGGCAGGCAGGCAGG".to_string(),
            "CAGGCAGGCAGGCAGGCAGGCAGGCAGGCAGACAGGCAGCCAGGCAGGCAGGCAGG".to_string(),
            "CAGGCAGGCAGGCAGGCAGGCAGGCAGGCAGACAGGCAGCCAGGCAGGCAGGCAGG".to_string(),
            "CAGGCAGGCAGGCAGGCAGGCAGGCAGGCAGACAGGCAGCCAGGCAGGCAGGCAGG".to_string(),
            "CAGACAGGCAGGCAGGCAGGCAGGCAGGCAGACAGGCAGCCAGGCAGGCAGGCAGG".to_string(),
            "ACAGACAGGCAGGCAGGCAGGCAGGCAGGCAGACAGGCAGCCAGGCAGGCAGGCAGACAGAA".to_string(),
            "CAGGCAGGCAGGCAGGCAGGCAGGCAGGCAGACAGGCAGCCAGGCAGGCAGGCAGG".to_string(),
            "CAGGCAGGCAGGCAGGCAGGCAGGCAGGCAGACAGGCAGCCAGGCAGGCAGGCAGG".to_string(),
            "CAGGCAGGCAGGCAGGCAGGCAGGCAGGCAGACAGGCAGCCAGGCAGGCAGGCAGG".to_string(),
            "CAGGCAGGCAGGCAGGCAGGCAGGCAGGCAGACAGGCAGCCAGGCAGGCAGGCAGG".to_string(),
        ];
        let cons = consensus(
            &seqs,
            10,
            20,
            &crate::repeats::RepeatInterval {
                chrom: "chr1".to_string(),
                start: 1,
                end: 100,
                created: None,
            },
            PoaScoring::default(),
        );
        println!("Consensus: {}", cons.seq.unwrap());
        println!("Num reads: {}", cons.support);
        println!("Std dev: {}", cons.std_dev);
        println!("Consensus score: {}", cons.score);
    }

    #[test]
    fn test_consensus_2() {
        let seqs = vec![
            "GGGGGGAGGAGGGGGGAGGAGGGGGGAGGAGGACAGGCAAGGAGGCCTCGAAGAGGGAAGACCAGAGAAGGAGGGGAAGGGGAGAGCGGAGGGAAAGAGGTGGGATAAGGAAGAGAAGGAGGAGGAAGGGGAAGAGGGAGG".to_string(),
            "GGGGGGAGGAGGGGGGAGGAGGGGGGAGGAGGACAGGCAAGGAGGCCTCGAAGAGGGAAGACCAGAGAAGGAGGGGAAGGGGAGAGCGGAGGGAAAGAGGTGGGATAAGGAAGAGAAGGAGGAGGAAGGGGAAGAGGGAGG".to_string(),
            "AGGAGGAGGGAGGCGGGGGGAGGAGGGGGGAGGAGGGGGGAGGAAGGACAGGCAAGGAGGCCTCGAAGAGGGAAGACCAGAGAAGGAGGGGAAGGGGAGAGCGAGAAGAAGGTGGGATGAGGAAGGAAAGGAGC".to_string(),
            "GGGGGAGGAGGGGGGAGGAGGGGGGAGGAGGACAGGCAGTGGTGGCTCCGAAGAGGGAAGACCAGAGAAGGAGGGGAAGGGGAGAGCGGAGGGAAAGAGGTGGGATAAGGAAGAGAAGGAGGAGGAAGGGGAAGAGGGAGG".to_string(),
            "GGGGGGAGGAGGGGGGAGGAGGGGGGAGGAGGGACAGGGCAAGGAGGCCTCGAAGAGGGAAGACCAGAGAAGGAGGGGAAGGGGAGAGCGAGGGAGAGAAGGATGAGGGAGTAGGGAAGGAAAGGAGGAGGAAGGGGAAGAGGGAGG".to_string(),
            "GGGGGGAGGAGGGGGGAGGAGGGGGAGGAGGACAGGCAAGGAGGCCTCGAAGAGGGAAGACCAGAGAAGGAGGGGAAGGGGAGAGCGGAGGGAAAGAGGTGGGATAAGGAAGAGAAGGAGGAGGAAGGGAAGAGGGAGG".to_string(),
            "GGGGGGAGGAGGGGGAGGAGGGGAGGAGGACAGGCAAGGAGGCCTCGAAGAGGGAAGACCAGAGAAGGAGGGGAAGGGGAGAGCGGAGGGAAAGAGGTGGGATAAGGAAGAGAAGGAGGAGGAAGGGGAAGAGGGAGG".to_string(),
            "GGGGGGAGGAGGGGGGAGGAGGGGGGAGGAGGACAGGCAAGGAGGCCTCGAAGAGGGAAGACCAGAGAAGGAGGGGAAGGGGAGAGCGGAGGGAAAGAGGTGGGATAAGGAAGAGAAGGAGGAGGAAGGGGAAGAGGGAGG".to_string(),
            "GGGGGGAGGAGGGGGAGGAGGGGGGAGGAGGACAGGCAAGGAGGCCTGGGAGAGAGGGAAGACCAGAGAAGGAGGGGAAGGGGAGAGCGAGGAAAGAGGTGGGATAAGGAAGAAGGAGGAGGAAGGGGAAGAGGGAGG".to_string(),
            "GGGGGGAGGAGGGGGGAGGAGGGGGAGGAGGACAGGCAAGGAGGCCTCGAAGAGGGAAGACCAGAGAAGGAGGGGAAGGGGAGAGCGGAGGGAAAGAGGTGGGATAAGGAAGAGAAGGAGGAGGAAGGGGAAGAGGGAGG".to_string(),
            "GGGGGGAGGAGGGGGGAGGAGGGGGAGGAGGACAGGCAAGGAGGCCTCGAAGAGGGAAGACCAGAGAAGGAGGGGAAAGGGGAGAGCGGAGGGAAGAGGTGGGATAAGGAAGAGAAGGAGGAGGAAGGGGAAGAGGGAGG".to_string(),
            "GGGGGGAGGAGGGGGAGGAGGGGGGAGGAGGACAGGCAAGGAGGCCTCGAAGAGGGAAGACCAGAGAAGGAGGGGAAGGGGAGAGCGGAGGGAAAGAGGTGGGATAAGGAAGAGAAGGAGGAGGAAGGGGAAGAGGGAGG".to_string(),
            "GGGGGGAGGAGGGGGAGGAGGGGAGGAGGACAGGCAAGGAGGCCTCGAAGAGGGAAGACCAGAGAAGGAGGGGAAGGGGAGAGCGGAGGGAAAGAGGTGGGATAAGAAGAGAAGGAGGAGGAAGGGGAAGAGGGAGG".to_string(),
            "GGGGGGAGGAGGGGGGAGGAGGGGGGAGGAGGACAGGCAAGGAGGCCTCGAAGAGGGAAGACCAGAGAAGGAGGGGAAGGGGAGAGCGGAGGGAAAGAGGTGGGATAAGGAAGAGAAGGAGGAGGAAGGGGAAGAGGGAGG".to_string(),
            "GGGGGAGGAGGGGGGAGGAGGGGGGAGGAGGACAGGCAAGGAGGCCTCGAAGAGGGAAGACCAGAGAAGGAGGGGAAGGGGAGAGCGGAGGGAAAGAGGTGGGATAAGGAAGAGAAGGAGGAGGAAGGGGAAGAGGGAGG".to_string(),
            "GGGGGGAGGAGGGGGAGGAGGGGGAGGAGGACAGGGAAGGAGGCCTCGAAGAGGGAAGACCAGAGAAGGAGGGGGAAGGGAGAGCGGAGGGGAAAAGAGGTGGGATAAGGAAGAGAAGGAGGAGGAAGGGGAAGAGGGAGG".to_string(),
        ];
        let cons = consensus(&seqs, 1, seqs.len(), &dummy_repeat(), PoaScoring::default());
        assert!(cons.seq.as_ref().is_some_and(|s| !s.is_empty()));
        println!("Consensus: {}", cons.seq.unwrap());
        println!("Consensus score: {}", cons.score);
    }

    #[test]
    fn test_consensus_3() {
        let seqs = vec![
            "TCTTTCTTTCTTTCTTTCCTTTCCTTTCCTTTCCTTTCCTTCCTTTCCTTCCTTCCTTCCTTCCTTCCTTCCTCCCTCACTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCCCTCTCCCTCTCCCTCTCCCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCCCTCCCTCTCCCTCTCTCTCTCTCTCTCCCTCTCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCTCTCCTCTCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCTCTCCTCTCCTCTCTCTCCTCTCCTCTCCTCTCCCTCTCCCTCTCCTCTCCCTCTCCCTCTCCTGCTCCCTCTCCTCTCCCTCTCCCTCCTCCTCCCTCTCCTCTCCCTCCCTCCTTTCTCCCTCTCTCCCTCTCCCTCTCCCTCTCCCCTCTCCTCTCCTCTCCCTCTCCCTCTCCCCTCTCCCTCTCCCTCCTCCCTCTCCTCTCCCTCTCCCTCTCCTCCTCCCTCTCCCTCTCCCTCTCTCCCTCTCCCTCTCCCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCCTCTCCCTCTCCCCTCTCCTCTCCTCTCCCCCTCTCCTCTCCTCTCCCTCTCCCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTTTCTTTCTTT".to_string(),
            "TTCTTTCCTTCCTTTCCTTTCTTTCCTTTCCTTCCTTCCTTTCCTTCCTTCCTTCCTTCCTTCCTTCCTCCCTCCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCCCTCTCCCTCTCCCTCTCCCTCCCTCTCCCTCTCTCCCTCTCCCTCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCTCTCTCTCTCTCTCTCTCTCCCTCTCCCTCCCTCTCTCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCTCTCTCTCTCTCTCCCTCCTCTCCCTCTCCCTCTCCCTCTCTCCCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCCCTCTCCCTCTCCCTCCCTCTCCCTCCCTCTCCCTCTCCCTCTCCCTCTCTCCCTCTCCCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCCCTCCCTCCCTCCCTCTCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCTCCCTCCTCCCTCTCTCTCTCTCCCTCTCCCTCCCTCTCCCTCTCCCTCTCCCTCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTTTCTTTCTTTCTTTCTTCCC".to_string(),
            "CTTTTCTTTCCTTTCCTTTCCTTTCCTTTCCTTTCCTTCCTTTCCTTCCTTCCTTCCTTCCTTCCTTCCCTCCCTCCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCCTTTCTCCTCTCCCTCCCTCCCTCTCCCTCTCCCCCCTCCCTCTCCCTCCCTCTCCCTCTCTCCCTCTCTCTCCCTCCTCCTCTCCTCCCTCCCTCCTCTCCTCTCCCTCCTCCTCTCCTCTCCTCTCTCCTCCTCTCCCTCTCCTCTCCTCCCTCCTCTCCTCCTCCCTCTCCTCTCCTCCTCTCCTCTCCTCTCCTCTCCTCCTCTCCTCTCTCCCTCCTCTCCTCCCTCCTCCCTCCCTCCTCCTCTTCCCTCCTCTCCTCCTCCTCCTCCTCTCCTCTCCTCCTCCTCCCTCCCTCCTCTTCCTCTCTCCCTCTCCTTCTCCCTCTCTCTCCCTCCCTCCTCTCCTCCTCTCCCTCTCCTCCTCTCCCTCTCCTCTCCTCCTCCTCTCCTCTCCTCTCCTCTCCTCTCCTCTCCCTCTCCCTCTCCCTCTCCCTCCCTCCCTCTCCCTCTCCCTCTCTCTCCCTCCCTCTCCCTCCCTCTCCTCTCTCTCTCTCTCTTCTTTTTCTTTCTTTCCTCG".to_string(),
            "TCTTTCTTTCTTTCTTTTCCTTTCCTTTCCTTTCCTTTCCTTTCCTTCCTTTCCTTCCTTCCTTCCTTCCTTCCTTCCTCCCTCCCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCCCTCTCCCTCTCCCTCTCCCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCCCTCTCCCTCTCCCTCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCTCTCCCTCTCCCTCTCTCCCTCTCCCTCTCCCTCTCCCTCTCTCTCTCCCTCTCCCTCTCCCTCTCCCTCTCTCTCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCTCTCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCTCTCTCCCTCCCTCCCTCTCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCTCCCTCCCTCTCTCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTTTCTTTCTTT".to_string(),
            "TCTTTCTTTCCTTTCCTTTCCTTTCTTTCTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTCCCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCCCTCTCCCTCTCCCTCTCCCTCCCTCTCCCTCCCTCTCCCTCTCCCTCTCCCTCTCTCCCTCTCCCTCTCCCTCTCCCTCCTCTCCCTCTCCTCTCCCTCTCTCTCTCCCTCTCCCTCTCCCTCCCTTCCTCCCTCTCCTCCTCCCTCTCTCTCCTCTCCCTCTCTCTCTCTCCCTCTCCCTCTCCCTCTCTCTCTCTCCCCTCTCCCTCTCCCTCCTCCCTCTCCCTCTCCTTCTCTCTCCCTCTCCTCTCCCTCTCCCTCTCCCCTCCCTCTCCCTCTCCCTCTCCCTTCTCCCTCCTCTTCTCCTCTCCTCTCCTCTCCTCTCCCTCCTCCTCTCCTCCTCTCCTCTCCCTCTCCTCTCCCTCTCCCTCTCCTCTCCCCTCTCCTCTCCCTCTCCCTCTCCCTCTCCTCTCCCTCTCCCTCTCCTCTCCCTCTCCTCTCCTCTCCCTCTCCCTCTCCCTCCTCCCTCTCCCTCTCCTCCCTCTCCCTCCTCCCTCTCCCTCTCCTCCTCCCTCCCTCCCTCCCTCTCCCCTCTCCCTCTCCCTCCCTCCCTCTCCTCTCCCTCTCCTCTCCTCTCCTCCTCCCTCTCCCTCCCTCCCTCTCTCCCTCTCTCTCTCTCTCTCTCTCTTTCTTTCTTTCTTTTTC".to_string(),
            "TCTTTCTTTCTTTCTTTCCTTTCCTTTCCTTTCTTTCCTTTCCTTCCTTTCCTTCCTTCCTTCCTTCCTTCCTTCCTCCCTCCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCCCTCTCCCTCTCCCTCTCCCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCCCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCTCTCCCTCTCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCCTCTCCTCTCCCTCTCCCTCTCCCTCTCCTCTCCCCTCTCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCTCTCCTCTCCTCTCCTCTCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTTTCTTTCTTT".to_string(),
            "TCTTTCTTTCTTTCTTTCCTTTCCTTTCCTTTCCTTTCCTTTCCTTCCTTTTCCTTCCTTCCTTCCTTCCTTCCTTCCTCCCTCCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCCCTCTCCCTCTCCTCTCCCTCCCTCTCCCTCTCCCTCTCCCTCCTCCCTCTCCCTCTCCCTCTCCCTCTCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCCCTCTCCCTCTCCCTCTCCCTCTCCTCTCTCTCCCTCTCCCTCTCCCTCTTTCCCTCTCCCTCTCCCTCTCCCCTCTCCCTCTCCCTCTCCCTCCCTCTCCCTCCCTCTCCCTCTCCCTCCTCTCTCTCCCTCTCCCTCCTCCCTCTCCTCCTCCTCTCCCTTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTTCTCCCTCTCCCTCCTCTCCCTCTCCCTCTCTCTCCCTCTCCTCTCCCTCCCTCCCTTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCTCTCCCTCCCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCTCTCCTCTCCCTCTCCCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTTTCTTTCTTT".to_string(),
            "TCTTCCTTTTCCTTCTTTCCTTTTCCTTTTCCTTCCTTCCTTCCTCCTTCCTTCCTTCCTTCCTCCCTCCCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCCCTCCCTCCCTCCCTCTCCCTCCCTTCTCCCTCTCCCTTCTTCCTTCCTTTTCTTCCTTTTTCTTCCTTCTTTCCTTCTTCCTTTCTTCTTTCTTTCCTTTCCTTTCTTTCTTTCTTTCTTCTTTTTCTTCTTTCTTTCTTTTTTCTTTCTTTCTTCCTTTTTTCTTTCTTTTTCTTTCCTTTTTCTTCTTTTTCTTTCTTTCTTTCTTCTTTTCTTTCTTTCTTTCTTTCTTTCCTTTTTCTTCCTTTCTTTCTTTCCTTTTTCTTTCTTTCTTCCTTCTTTCTTTCTTTCTTTCTTTTTCTTTCCTTCCTTTTCTTTCTTTCTTTCTTTCTTTCTTTCCTTTTTCTTCTTTCCTTCTCTCTTTTCTTCCTTTTTCCTTTCTTTCCTTCTTTTTCTTTCCTTTCTTTTTCTTCCTTTCTCTCTTTTTCTTTCTTTCCTTTTCTTTCCTTTCTCTCTTTCTTCCTTTTTCTTCCTTTTCTTTCTTCTTCTTTCTTTCCTTTTCTTTCCTTTCTTTCCTTTTTCTTCCTTTTCTTTCTTTCTTTCTTTCTTCCTTCTTTCTTTCCTTTCTTTCTTTCTTCTTTCTTTCCTTTCTTTCTTTCCTTTTTCTTCTTTCCTTCTTCTTTCCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTTTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCCCTCTCCCTCTCCCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTTTCTTCTTTCTTTCTTTCCTG".to_string(),
            "TCTTTCTTTCTTTCTTTGCTTTCCTTTCTTTCCTTTCCTTTCCTTTCCTTCCTTCCTTCTTCCTTCCTTCCCTTCCTCCCTCCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCCCTCTCCCTCTCCCTCTCCCTCCCTCTCCCTCTCCCTCTCCCTCTCTCTCCCTCTCTCTGCCTCTCCCTCTCCCTCTCCCTCTCCTCTCCCTCTCCCTCCCTCCCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCTCTCCCTCCCTCCCTCTCTCTCCTCTCCCTCACCCTCTCCCTCTCCCTCTCCCTCTCCTCTCCCTCTCCCTCTCCTCTCCCTCCCTCTCCCCTCCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCTCTCCCTCTCCCTCTCCCTCTCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCTCTCCCTCTCCCTCTCCCCTCCTCCTCCCTCTCCTCTCCCTCTCCCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCTCTCCCTCTCCCTCTCCCTCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCCTCCCTCTCCCTCTCCCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTTTCTTTCTTT".to_string(),
            "TCTTTCTTTCTTTCTTTCCTTTCCTTTCCTTTCCTTTCCTTTCCTTCCTTTCCTTCCTTCCTTCCTTCCTTCCTTCCTCCCTCCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCCCCTCTCCCTCTCCTCCTCCTCCTCTCCTCTCCTCTCCCTCTCCTCTCCTCTCCCTCTCTCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCCTCTCCCTCTCCTCTCCCTCTCCTCTCCTCTCCCTCTCCTCCTCTCCCTCTCCTCTCCCTCTCCTCTCCCCTCTCCTCTCCCTCTCCCTCTCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCTCTCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCTCTCCCTCTCCCCTCTCCCTCTCCCTCTCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCTCTCCTCTCCTCTCCTCTCCCTCTCCTCTCCCTCTCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCCCTCTCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTTTCTTTCTTT".to_string(),
            "TCTTTCTTTTCTTTCTTTCCTTTCCTTTCCTTTCTTTCCTTTCCTTCCTTTCCTTCCTTCCTTCCTTCCTTCCTTCCTCCCTCCCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCCCTCTCCCTCTCCCTCCCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCTCCCTCTCCCTCTCTCCCTCTCCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCTCCCTCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCTCTCCCTCTCCCTCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCTCCCTCTCCTCTCTCCCTCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCTCCCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCTCTCTCTCTCCCTCTCCCTCTCCCTCTCCCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTTTCTTTTTCTTT".to_string(),
            "TCTTTCTTTCTTTCTGCTTTCCTTTCCTTTCCTTTCCTTTCCTTTCCTTCCTTTCCTTCCTTCCTTCCTTCCTTCCTTCCTCCCTCCTCTCTCTCTCTCTCTCTCTCTCTCTTTCCTCCTCCCTCCCTCCCTCCCTCTCCCTCCCTCTCCCTCCTCCCTCTCCCTCTCCCTCCTCTCCCTCCCTCCCTCTCCCTCCCTCTCCCTCTCCCTCTCTCCCTCTCCTCTCCTCTCCCTCTCCCTCTCCTCTCCCTCCTCCCTCTCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCTCTCCCTCTCTCCCACCTCTCCCTCTCCCTCTCCCTCCCTCCCTCTCCTCTCCCTCTCCCTCCCTCCTCCCTCTCCCTCTCCCTCTCCCTCTCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCTCTCCTCTCCCCCTCTCCCTCTCCCTCTCCCTCCTCCCTCTCCCCTCTCCCCCACCCCTCTCCCCTCTCCCTCCCTCTCCCTCTCCCCTCTCCCTCTCCCTCTCCCTACTCCCTCCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCTCCCTCTCCCTCTCCCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTTTCTTTCTTT".to_string(),
            "TCTTTCTTTCTTTCTTTCCTTCCTTTCCTTTTCCTTTCCTTTCCTTCCTTTCCTTCCTTCCTTCCTTCCTTCCTTCCTCCCTCCCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCCCTCTCCCTCCCTCTCTCCCTCCCTCTCCCTCTCTCCCTCCCTCTCCCTCTCCCTCTCCCTCTCTCTCTCCCTCCCTCTCCCTCTCCCTCTCTCTCCCTCCCTCTCCCTCCCTCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCTCCCTCCCTCTCCCTCTCCCTCTCTCCCTCTCCCTCTCCCTCTCCCTCTCTCCCTCTCCCTCCCTCTCCCTCTCCCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCTCCCTCTCCCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCCCTCTCCCTCTCTCTCTCTCTCTCCCTCTCCCTCTCTCCCTCCCTCCCTCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCCTCTCCCTCTCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTTTCTTTCTTT".to_string(),
            "TCTTTCTTTCTTTCTTTCCTTTCCTTTCCTTTCCTTTCCTTTCCTTCCTTTCCTTCCTTCCTTCCTTCCTTCCTCCCTCCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCCCTCTCCCTCTCCCTCTCCCTCCCTCTCCCTCTCCCTCCTCTCCCTCTCCCTCTCCCTCTCCCTCCTCTCCTCTCCTCTCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCTTCTCCTCTCCTCTCCCTCTCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCCCCTCTCCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCTCTCCCTCTCCCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCTCCCTCTCCCTCCTCTCCCTCTCCCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTTTCTTTCTTT".to_string(),
            "CTTTCTTTCTTTCTTTCTTTCCTTTCCTTTCCTTTCCTTTCCTTCCTTTCCTTCCTTCCTTCCTTCCTTCCTTCCTCCCTCCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCCCTCTCCCTCTCCCTCTCCCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCCCTCCTCTCCCTCCTCCTCTCCCTCTCCCTCTCCCTCTCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCTCTCCCTCTCCCTCTCCCTCCTCCCTCCCTCCCTCTCCCTCTCCCTCCTCCCTCTCCTCTCCTCTCCCTCTCCCTCTTTCCCTCTCCTCTCCCTCTCCTCTCCCTCTCCCCTCTCCCTCTCCTCTCCTCTCCCTCTCCTCTCCTCTCCCTCTCCTCTCCCTCTCCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCTCCCTCCCTCTCCTCTCCCTCTCCCTCTCCTCTCCCTCTCCTCTCCCTCTCCCTCTCCCCTCTCCCTCTCCCCTCTCCCTCTCCTCTCCTCTCCCTCTCCCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTTTCTTTCTTTTCTTTC".to_string(),
                ];
        let cons = consensus(&seqs, 1, seqs.len(), &dummy_repeat(), PoaScoring::default());
        assert!(cons.seq.as_ref().is_some_and(|s| !s.is_empty()));
        println!("Consensus: {}", cons.seq.unwrap());
        println!("Consensus score: {}", cons.score);
    }

    #[test]
    fn test_consensus_4() {
        let seqs = vec![
            "TTTCTTTCTTTCTTTCTTTCTT".to_string(),
            "TTTCTTTCTTTCTTTCTTTCTTTCTTTTTTTT".to_string(),
            "TTTCTTTCTTTCTTTCTTTCTTTCTTTTT".to_string(),
            "TTTCTTTCTTTCTTTCTTTCTT".to_string(),
            "TTTCTTTCTTTCTTTCTTTCTTTCTTT".to_string(),
            "ACTTCTTTCTTTCTTTCTTTCTTTCTTTT".to_string(),
            "TTTCTTTCTTTTCTTTCTTTCTTTCTT".to_string(),
            "TTTTCTTTCTTTCTTTCTTTCTTTCTTTT".to_string(),
            "TTTCTTTTCTTTCTTTCTTTCTTTCTTTTT".to_string(),
            "TTTCTTTCTTTTCTTTCTTTCTTTCTTT".to_string(),
            "TTTCTTTCTTTCTTTCTTTCTTTCT".to_string(),
            "TTTCTTTCTTTCTTTCTTTCTTTCTTTT".to_string(),
            "CTCTCTTTCTTTTCTTTCTTTCTTTTTCTTT".to_string(),
            "TTTCTTTCTTTCTTTCTTTCTT".to_string(),
            "TTTCTTTCTTTCTTTCTTTTCTTTCTTTT".to_string(),
            "TCGTTTCTTTCTTTCTTTCTTTCTTTCTTT".to_string(),
            "TTCTTTCTTTCTTTCTTTTTTTTTC".to_string(),
            "TTTCTTTCTTTCTTTCTTTCTTTCTT".to_string(),
            "TTTCTTTCTTTTTTTTCTTTCTTTCTTTCTTTT".to_string(),
            "TTTCTTTCTTTCTTTCTTTCTTTCTTTTT".to_string(),
        ];
        let cons = consensus(&seqs, 1, seqs.len(), &dummy_repeat(), PoaScoring::default());
        assert!(cons.seq.as_ref().is_some_and(|s| !s.is_empty()));
        println!("Consensus: {}", cons.seq.unwrap());
        println!("Consensus score: {}", cons.score);
    }

    #[test]
    fn test_consensus_5() {
        let seqs = vec![
            "TATATATATATAAACATATATTATATATATAAAATATAACATATATAAACATATATATTATATATATA".to_string(),
            "TATATATATATAAACATATATTATATATGTAATATAAACATATATAAACATATATTATATATA".to_string(),
            "TATATATATATAAACATATATTATATATAATATAAACATATATAAACATATATTATATATATA".to_string(),
            "TATATATATATAAACATATATTATATATGTAATATAAACATATATAAACATATATTATATATATA".to_string(),
            "ATATATATAAACATATATTATATATGTAATATAAATATATATAAACATATATTTATATATATA".to_string(),
            "TATATATATATAAACATATATTATATATGTAATATAAACATATATGTATACATATATATACA".to_string(),
            "TATATATATATAAACATATATTATATATGTAATATAAATATATATAAACATATATTATATATATA".to_string(),
            "TATATATATATAAACATATATTATATATGTAATATAAACATATATAAACATATATTATATATATA".to_string(),
            "TATATATTTATAAACATATATTATGTATGTAATATAAACATATATAAACATATATTATATATA".to_string(),
            "TATATATATATAAACATATATTATATATATATAAACATATATAAACATATATTATATATATA".to_string(),
            "TATATATATATAAACATATATTCTATATATGTAATATAAACATATATAAACATATATTATCTATATA".to_string(),
            "TATATATATATAAACATATATTATATATAATATAAACATATAAACATATATTATATATATA".to_string(),
            "TATATATATATAAACATATATTATATATGTAATATAAACATATATAAACATATATTATATATATA".to_string(),
            "TATATATATATAAACATATATTATATATGTAATATGTTTTCTATATGTTGCTATATTATACAACATA".to_string(),
            "ATATATATATATAAACATATATTATATATGTAATATAAACATATATAAACATATATTATATATATATA".to_string(),
            "TATATATATATAAACATATATTATATATGTAATATAACATATATAAACATATATTATATATATA".to_string(),
        ];
        let cons = consensus(&seqs, 1, seqs.len(), &dummy_repeat(), PoaScoring::default());
        assert!(cons.seq.as_ref().is_some_and(|s| !s.is_empty()));
        println!("Consensus: {}", cons.seq.unwrap());
        println!("Consensus score: {}", cons.score);
    }
}
