use std::collections::{BTreeMap, HashMap};

// write a function that takes a long string like CAGCAGCAGCAGCGGCGGCGGCAGCAGCAG and converts it to a condensed representation like (CAG)4(CGG)3(CAG)3
// the function should identify repeated parts and replace them with (repeat)number

#[allow(dead_code)]
fn create_motif(_seq: &str) -> String {
    unimplemented!()
}

/// Longest period `dominant_motif` considers: homopolymers up to VNTR units such as TYMS's
/// 28-mer or MUC1's 60-mer. A unit longer than the longest period considered would be
/// misreported as whichever shorter period happens to match, so this is generous; the cost
/// is bounded by half the allele length, since a unit must occur twice.
pub const MOTIF_MAX_PERIOD: usize = 100;
/// A stretch where the sequence matches itself shifted by the period scores +1 per matching
/// base and minus this per mismatch, so it extends through sequence that is at least 67%
/// pure and ends at flanking sequence. Catalog intervals (adotto's in particular) pad the
/// repeat with flank and can hold several repeats, so the whole interval cannot be scored.
/// 67% rather than 75% keeps codon repeats such as polyalanine's GCN, whose third base
/// varies, in one stretch; it also found a unit in more of adotto's chr2 intervals (75% vs
/// 72%) while agreeing at least as often with their annotated motifs.
const MOTIF_MISMATCH_PENALTY: i32 = 2;
/// Minimum score of the best stretch for an allele to count as a tandem repeat at all:
/// an 11-base homopolymer, (CA)6 or (CAG)4 and a bit. Random sequence does not get there.
const MOTIF_MIN_SCORE: i32 = 10;
/// Multiples of the true period score nearly as well as the period itself, and better
/// across an interruption: the smallest period scoring at least this fraction of the best
/// is the unit, so ((CAG)5 CAA)n is CAG rather than the 18-mer.
const MOTIF_PERIOD_MARGIN: f64 = 0.85;
/// An interrupted or impure repeat can still score best at a multiple of its unit, which
/// spans the irregularity. A divisor of the chosen period whose stretch scores at least this
/// fraction of the chosen one (and `MOTIF_MIN_SCORE`) is preferred: (GCC)n with interruptions
/// is GCC. On the STRchive catalog this turned two composite units back into codons and
/// did not cost agreement with adotto's annotations.
const MOTIF_DIVISOR_MARGIN: f64 = 0.5;

/// The dominant repeat unit of an allele sequence, or None when it holds no tandem repeat of
/// a unit up to `MOTIF_MAX_PERIOD` (see `MOTIF_MIN_SCORE`).
///
/// Per period, the best-scoring stretch where the sequence matches itself shifted by that
/// period is found (see `MOTIF_MISMATCH_PENALTY`); the period is the smallest one whose
/// stretch scores close to the best (see `MOTIF_PERIOD_MARGIN`), or a divisor of it that
/// scores reasonably (see `MOTIF_DIVISOR_MARGIN`). Codon repeats with a degenerate base,
/// such as polyalanine's GCN, are reported as their most common codon. The unit is the most
/// frequent `period`-mer within that stretch, with all rotations counted together (CAG, AGC
/// and GCA are one repeat read in different phases). It is reported in the rotation most
/// common in the reference's frame, i.e. in period-sized chunks from the start of the
/// reference's own stretch, when the reference holds the same unit, else in the allele's.
/// Anchoring on the reference makes every allele and sample at a locus report the same
/// rotation: HTT gives the same unit for both alleles, rather than CAG for one and AGC for
/// the other. Only the dominant unit is reported: an interrupted or mixed allele such as
/// (AAAAG)n(AAGGG)m gives whichever unit covers more of it.
pub fn dominant_motif(seq: &str, reference: &str) -> Option<String> {
    let seq = seq.as_bytes().to_ascii_uppercase();
    let (period, unit, stretch) = repeat_unit(&seq)?;
    let reference = reference.as_bytes().to_ascii_uppercase();
    let reference_stretch = match reference.len() > period {
        true => best_stretch(&reference, period).1,
        false => 0..0,
    };
    let rotation = rotation_in_frame(&reference[reference_stretch], period, &unit)
        .or_else(|| rotation_in_frame(&seq[stretch], period, &unit))
        .unwrap_or(unit);
    Some(String::from_utf8_lossy(&rotation).into_owned())
}

/// (period, unit as its least rotation, the stretch it was found in), see `dominant_motif`.
fn repeat_unit(seq: &[u8]) -> Option<(usize, Vec<u8>, std::ops::Range<usize>)> {
    let max_period = MOTIF_MAX_PERIOD.min(seq.len() / 2);
    let stretches: Vec<(i32, std::ops::Range<usize>)> =
        (1..=max_period).map(|p| best_stretch(seq, p)).collect();
    let best = stretches.iter().map(|(score, _)| *score).max()?;
    if best < MOTIF_MIN_SCORE {
        return None;
    }
    let index = stretches
        .iter()
        .position(|(score, _)| *score as f64 >= best as f64 * MOTIF_PERIOD_MARGIN)?;
    let candidate = index + 1;
    let candidate_score = stretches[index].0 as f64;
    let period = (1..candidate)
        .filter(|d| candidate % d == 0)
        .find(|d| {
            let score = stretches[d - 1].0;
            score >= MOTIF_MIN_SCORE && score as f64 >= candidate_score * MOTIF_DIVISOR_MARGIN
        })
        .unwrap_or(candidate);
    let stretch = stretches[period - 1].1.clone();

    // count raw k-mers first: far fewer distinct ones than windows in a repeat, so each
    // distinct k-mer is rotated once rather than once per position
    let mut raw: HashMap<&[u8], usize> = HashMap::new();
    for window in seq[stretch.clone()].windows(period) {
        *raw.entry(window).or_default() += 1;
    }
    let mut classes: BTreeMap<Vec<u8>, usize> = BTreeMap::new();
    for (kmer, count) in raw {
        *classes.entry(least_rotation(kmer)).or_default() += count;
    }
    // the BTreeMap order makes ties deterministic: the lexicographically smallest unit wins
    let (unit, _) = classes.into_iter().rev().max_by_key(|(_, count)| *count)?;
    Some((period, unit, stretch))
}

/// The best-scoring stretch where `seq` matches itself shifted by `period` (maximum
/// subarray, Kadane), as its score and the bases it covers, both copies included.
fn best_stretch(seq: &[u8], period: usize) -> (i32, std::ops::Range<usize>) {
    let (mut best, mut best_range) = (0, 0..0);
    let (mut current, mut start) = (0, 0);
    for i in 0..seq.len() - period {
        if current <= 0 {
            current = 0;
            start = i;
        }
        current += if seq[i] == seq[i + period] {
            1
        } else {
            -MOTIF_MISMATCH_PENALTY
        };
        if current > best {
            best = current;
            best_range = start..i + 1 + period;
        }
    }
    (best, best_range)
}

/// The rotation of `unit` seen most often among the period-sized chunks of `seq` taken from
/// its start, or None when no chunk is a rotation of `unit`.
fn rotation_in_frame(seq: &[u8], period: usize, unit: &[u8]) -> Option<Vec<u8>> {
    let mut in_frame: BTreeMap<&[u8], usize> = BTreeMap::new();
    for chunk in seq
        .chunks_exact(period)
        .filter(|c| least_rotation(c) == unit)
    {
        *in_frame.entry(chunk).or_default() += 1;
    }
    let (rotation, _) = in_frame.into_iter().rev().max_by_key(|(_, count)| *count)?;
    Some(rotation.to_vec())
}

/// The lexicographically smallest rotation, in linear time (two-pointer minimum expression),
/// so that units up to `MOTIF_MAX_PERIOD` stay cheap to compare.
fn least_rotation(s: &[u8]) -> Vec<u8> {
    let n = s.len();
    let (mut i, mut j, mut k) = (0, 1, 0);
    while i < n && j < n && k < n {
        let (a, b) = (s[(i + k) % n], s[(j + k) % n]);
        if a == b {
            k += 1;
            continue;
        }
        if a > b {
            i += k + 1;
        } else {
            j += k + 1;
        }
        if i == j {
            j += 1;
        }
        k = 0;
    }
    let start = i.min(j);
    [&s[start..], &s[..start]].concat()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dominant_motif_period_and_phase() {
        let m = |seq: &str| dominant_motif(seq, "");
        // without a reference to anchor on, the allele's own frame
        assert_eq!(m(&"CAG".repeat(20)).as_deref(), Some("CAG"));
        assert_eq!(m(&"AAGGG".repeat(30)).as_deref(), Some("AAGGG"));
        assert_eq!(m(&"A".repeat(15)).as_deref(), Some("A"));
        // a homopolymer-rich unit is not mistaken for its homopolymer
        assert_eq!(m(&"AAAAAG".repeat(10)).as_deref(), Some("AAAAAG"));
        assert_eq!(m(&"GGCGCGGAGC".repeat(5)).as_deref(), Some("GGCGCGGAGC"));
        // a VNTR unit is found whole, not as a shorter period inside it
        let tyms = "CCGCGCCACTTGGCCTGCCTCCGTCCCG";
        assert_eq!(m(&tyms.repeat(3)).as_deref(), Some(tyms));
        // interruptions; the larger part of a mixed allele wins, in frame
        let interrupted = format!("AG{}CAA{}", "CAG".repeat(12), "CAG".repeat(8));
        assert_eq!(m(&interrupted).as_deref(), Some("AGC"));
        let mixed = format!("{}{}", "AAAAG".repeat(5), "AAGGG".repeat(40));
        assert_eq!(m(&mixed).as_deref(), Some("AAGGG"));
        assert_eq!(m("cagcagcagcagcag").as_deref(), Some("CAG"));
        // an interrupted allele is CAG, not the multiple that spans the interruption; even a
        // CAA in every third codon reads as its most common codon, as GCN repeats do
        assert_eq!(m(&format!("{}CAA", "CAG".repeat(5)).repeat(3)).as_deref(), Some("CAG"));
        assert_eq!(m(&"CAGCAGCAA".repeat(4)).as_deref(), Some("CAG"));
        assert_eq!(m(&"GGGGCC".repeat(8)).as_deref(), Some("GGGGCC"));
        // the repeat inside flanking sequence, as in a padded catalog interval
        let flanked = format!("TACAGAAGTTTTTCTAGTATTGC{}ATAAGTTAAATAAAATTAAAT", "CAAA".repeat(8));
        assert_eq!(m(&flanked).as_deref(), Some("CAAA"));
    }

    #[test]
    fn test_dominant_motif_takes_the_phase_of_the_reference() {
        let reference = "CAG".repeat(19);
        // an allele that starts out of phase still reports the reference's rotation
        assert_eq!(
            dominant_motif(&format!("AG{}", "CAG".repeat(30)), &reference).as_deref(),
            Some("CAG")
        );
        // RFC1: an expanded AAAAG allele read in another phase is still AAAAG
        let rfc1 = "AAAAG".repeat(11);
        assert_eq!(
            dominant_motif(&format!("AAG{}", "AAAAG".repeat(100)), &rfc1).as_deref(),
            Some("AAAAG")
        );
        // a unit absent from the reference keeps the allele's own frame
        assert_eq!(dominant_motif(&"AAGGG".repeat(100), &rfc1).as_deref(), Some("AAGGG"));
    }

    #[test]
    fn test_least_rotation() {
        assert_eq!(least_rotation(b"CAG"), b"AGC");
        assert_eq!(least_rotation(b"GAAAA"), b"AAAAG");
        assert_eq!(least_rotation(b"ABAB"), b"ABAB");
        assert_eq!(least_rotation(b"A"), b"A");
        // agrees with the quadratic definition
        for kmer in [&b"GGCGCGGAGC"[..], b"TTATG", b"CGCGCA", b"AACAAT"] {
            let naive = (0..kmer.len())
                .map(|r| [&kmer[r..], &kmer[..r]].concat())
                .min()
                .unwrap();
            assert_eq!(least_rotation(kmer), naive);
        }
    }

    #[test]
    fn test_dominant_motif_none_when_not_a_repeat() {
        assert_eq!(dominant_motif("", ""), None);
        assert_eq!(dominant_motif(".", "."), None);
        assert_eq!(dominant_motif("ACGTTGCAAGCTTAGCCATGAGTCAGTACGATCGAATCGGCTA", ""), None);
    }

    #[test]
    #[ignore]
    fn test_create_motif() {
        assert_eq!(create_motif("CAGCAGCAGCAGCGGCGGCGGCAGCAGCAG"), "(CAG)4(CGG)3(CAG)3");
    }
}
