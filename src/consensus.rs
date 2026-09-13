use bio::alignment::{
    pairwise::Scoring,
    poa::{Aligner, POAGraph},
};
use log::debug;
use petgraph::Incoming;
use petgraph::visit::Topo;
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

/// Scoring for the POA that builds an allele consensus.
///
/// rust-bio's POA never reads `Scoring::gap_extend` and charges `gap_open` for every gap
/// base, so the gap model is linear (<https://github.com/rust-bio/rust-bio/issues/677>).
///
/// These were once tunable, on the theory that a gap penalty cheap relative to a match let
/// a single read's insertion open its own node and be carried into the consensus. Measured
/// against the GIAB HG002 truth set that theory was half right: raising the gap penalty
/// from 12 to 70 did gain 11 points on its own, but *nothing* on top of `trim_fraction`
/// (-0.1), because both were addressing the same defect and trimming addresses it at the
/// source. The knobs were therefore retired rather than given new defaults; the values
/// below are the originals and should not be changed without re-measuring, since no
/// experiment has ever shown them to matter once the consensus endpoint is correct.
#[derive(Clone, Copy, Debug)]
pub struct PoaScoring {
    pub gap_open: i32,
    pub match_score: i32,
    pub mismatch: i32,
    /// Seed the POA graph with the read whose length is closest to the cluster median,
    /// rather than whichever read happened to be sampled first.
    pub medoid_seed: bool,
    /// Fraction of the cluster's reads an edge must carry for the consensus to run through
    /// it at either end. 0 disables trimming and uses rust-bio's own `consensus()`.
    pub trim_fraction: f64,
}

impl Default for PoaScoring {
    fn default() -> Self {
        Self {
            gap_open: -12,
            match_score: 3,
            mismatch: -4,
            medoid_seed: false,
            trim_fraction: 0.0,
        }
    }
}

/// Consensus through the POA graph, refusing to run out along poorly supported ends.
///
/// rust-bio's `Aligner::consensus()` picks its endpoint with `max_by_key` over a cumulative
/// score to which every edge contributes a weight of at least 1. The score therefore
/// increases strictly along every edge, so the argmax is *always a sink*: the consensus runs
/// to the deepest point in the graph rather than the best supported one. One read in fifteen
/// that extends a single base past the others creates a deeper sink and wins, with no vote
/// taken -- which is why the reported allele is one base too long on 8% of alleles at loci
/// that are truly homozygous reference, why the excess is 29:1 one-sided, and why it gets
/// *worse* with more reads instead of better.
///
/// The interior path choice upstream makes (heaviest edge first) is sound and is kept. Only
/// the ends are corrected: a terminal base is dropped unless the edge reaching it carries at
/// least `min_weight` reads. Reported upstream alongside rust-bio#677.
fn trimmed_consensus(graph: &POAGraph, min_weight: i32) -> Vec<u8> {
    let node_count = graph.node_count();
    // (weight of the best incoming edge, cumulative score, that predecessor's index)
    let mut best: Vec<(i32, i32, usize)> = vec![(0, 0, usize::MAX); node_count];
    let mut topo = Topo::new(graph);
    while let Some(node) = topo.next(graph) {
        let mut choice: (i32, i32, usize) = (0, 0, usize::MAX);
        for neighbour in graph.neighbors_directed(node, Incoming) {
            let index = neighbour.index();
            let weight: i32 = graph
                .edges_connecting(neighbour, node)
                .map(|e| *e.weight())
                .sum();
            let score = weight + best[index].1;
            if (weight, score, index) > choice {
                choice = (weight, score, index);
            }
        }
        best[node.index()] = choice;
    }

    let Some(end) = (0..node_count).max_by_key(|&i| best[i].1) else {
        return Vec::new();
    };
    let mut path = Vec::new();
    let mut pos = end;
    while pos != usize::MAX {
        path.push(pos);
        pos = best[pos].2;
    }
    path.reverse();

    // drop trailing nodes reached by an edge too few reads support, then leading nodes whose
    // edge into the rest is equally thin. The first node of the path has no incoming edge,
    // so it is judged by the edge leaving it.
    while path.len() > 1 && best[path[path.len() - 1]].0 < min_weight {
        path.pop();
    }
    let mut start = 0;
    while start + 1 < path.len() && best[path[start + 1]].0 < min_weight {
        start += 1;
    }

    path[start..]
        .iter()
        .map(|&i| graph.raw_nodes()[i].weight)
        .collect()
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

#[derive(Clone)]
pub struct Consensus {
    pub seq: Option<String>,
    pub support: usize,
    pub std_dev: usize,
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
    let seqs = remove_outliers(seqs, mean, std_dev, repeat);
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
        // I empirically determined the following parameters to be suitable,
        // but further testing on other repeats would be good
        // mainly have to make sure the consensus does not get longer than the individual insertions
        //
        // NOTE: rust-bio's POA aligner ignores the gap_extend term entirely (it never reads
        // `Scoring::gap_extend`) and applies a *linear* penalty of `gap_open` per gap base.
        // The -6 below is therefore inert; the effective gap model is -12 per base. Do not
        // bother tuning the second argument until/unless upstream POA gains affine gaps.
        // (reported upstream: https://github.com/rust-bio/rust-bio/issues/677)
        log::info!("Creating consensus for {repeat}");
        let (match_score, mismatch) = (poa.match_score, poa.mismatch);
        let scoring = Scoring::new(
            poa.gap_open,
            -6,
            move |a: u8, b: u8| {
                if a == b { match_score } else { mismatch }
            },
        );
        let seed = if poa.medoid_seed {
            medoid_index(&seqs_bytes)
        } else {
            0
        };
        debug!("{repeat}: seeding POA graph with read {seed} of {}", seqs_bytes.len());
        let mut aligner = Aligner::new(scoring, &seqs_bytes[seed]);
        for (i, seq) in seqs_bytes.iter().enumerate() {
            if i != seed {
                aligner.global(seq).add_to_graph();
            }
        }
        debug!("Added all sequences to graph");

        let consensus = if poa.trim_fraction > 0.0 {
            let min_weight = (poa.trim_fraction * seqs_bytes.len() as f64)
                .ceil()
                .max(1.0) as i32;
            trimmed_consensus(aligner.graph(), min_weight)
        } else {
            aligner.consensus()
        };
        debug!("Created consensus");
        let score = aligner.global(&consensus).alignment().score;
        debug!("Calculated score");
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
    repeat: &crate::repeats::RepeatInterval,
) -> Vec<&'a String> {
    // remove sequences that are shorter or longer than two standard deviations from the mean
    // except if the stdev is small
    debug!("{repeat}: mean: {}, std_dev: {}", mean, std_dev);
    if std_dev < 5 {
        debug!("std_dev < 5, not removing any outliers");
        seqs.iter().collect::<Vec<&String>>()
    } else {
        // avoid underflowing usize
        let min_val = mean.saturating_sub(2 * std_dev);
        let max_val = mean + 2 * std_dev;
        debug!("Removing outliers outside of [{},{}]", min_val, max_val);
        seqs.iter()
            .filter(|seq| seq.len() > min_val && seq.len() < max_val)
            .collect::<Vec<&String>>()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bio::alignment::pairwise::Scoring;

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
        let mut seqs_bytes = vec![];
        for seq in seqs.iter() {
            seqs_bytes.push(seq.to_string().bytes().collect::<Vec<u8>>());
        }

        // I empirically determined the following parameters to be suitable,
        // but further testing on other repeats would be good
        // mainly have to make sure the consensus does not get longer than the individual insertions
        let scoring = Scoring::new(-12, -6, |a: u8, b: u8| if a == b { 3 } else { -4 });
        let mut aligner = Aligner::new(scoring, &seqs_bytes[0]);
        for seq in seqs_bytes.iter().skip(1) {
            aligner.global(seq).add_to_graph();
        }

        let consensus = aligner.consensus();
        let score = aligner.global(&consensus).alignment().score;

        println!("Consensus: {}", std::str::from_utf8(&consensus).unwrap());
        println!("Consensus score: {}", score);
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
        let mut seqs_bytes = vec![];
        for seq in seqs.iter() {
            seqs_bytes.push(seq.to_string().bytes().collect::<Vec<u8>>());
        }

        // I empirically determined the following parameters to be suitable,
        // but further testing on other repeats would be good
        // mainly have to make sure the consensus does not get longer than the individual insertions
        let scoring = Scoring::new(-12, -6, |a: u8, b: u8| if a == b { 3 } else { -4 });
        let mut aligner = Aligner::new(scoring, &seqs_bytes[0]);
        for seq in seqs_bytes.iter().skip(1) {
            aligner.global(seq).add_to_graph();
        }

        let consensus = aligner.consensus();
        let score = aligner.global(&consensus).alignment().score;

        println!("Consensus: {}", std::str::from_utf8(&consensus).unwrap());
        println!("Consensus score: {}", score);
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
        let mut seqs_bytes = vec![];
        for seq in seqs.iter() {
            seqs_bytes.push(seq.to_string().bytes().collect::<Vec<u8>>());
        }

        // I empirically determined the following parameters to be suitable,
        // but further testing on other repeats would be good
        // mainly have to make sure the consensus does not get longer than the individual insertions
        let scoring = Scoring::new(-12, -6, |a: u8, b: u8| if a == b { 3 } else { -4 });
        let mut aligner = Aligner::new(scoring, &seqs_bytes[0]);
        for seq in seqs_bytes.iter().skip(1) {
            aligner.global(seq).add_to_graph();
        }

        let consensus = aligner.consensus();
        let score = aligner.global(&consensus).alignment().score;

        println!("Consensus: {}", std::str::from_utf8(&consensus).unwrap());
        println!("Consensus score: {}", score);
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
        let mut seqs_bytes = vec![];
        for seq in seqs.iter() {
            seqs_bytes.push(seq.to_string().bytes().collect::<Vec<u8>>());
        }

        // I empirically determined the following parameters to be suitable,
        // but further testing on other repeats would be good
        // mainly have to make sure the consensus does not get longer than the individual insertions
        let scoring = Scoring::new(-12, -6, |a: u8, b: u8| if a == b { 3 } else { -4 });
        let mut aligner = Aligner::new(scoring, &seqs_bytes[0]);
        for seq in seqs_bytes.iter().skip(1) {
            aligner.global(seq).add_to_graph();
        }

        let consensus = aligner.consensus();
        let score = aligner.global(&consensus).alignment().score;

        println!("Consensus: {}", std::str::from_utf8(&consensus).unwrap());
        println!("Consensus score: {}", score);
    }
}
