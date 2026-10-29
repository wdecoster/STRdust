use std::collections::{BTreeMap, HashMap};
use std::ops::Range;

// write a function that takes a long string like CAGCAGCAGCAGCGGCGGCGGCAGCAGCAG and converts it to a condensed representation like (CAG)4(CGG)3(CAG)3
// the function should identify repeated parts and replace them with (repeat)number

#[allow(dead_code)]
fn create_motif(_seq: &str) -> String {
    unimplemented!()
}

/// Longest period `allele_motifs` considers: homopolymers up to VNTR units such as TYMS's
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

/// Most units reported per allele.
const MOTIF_MAX_UNITS: usize = 3;
/// A further unit within a stretch must cover at least this share of it, and at least two
/// copies' worth of windows: a sequencing error only creates k-mers seen once or twice, so
/// a unit clearing both is a real part of the allele, such as AAAAG in (AAAAG)n(AAGGG)m.
const MOTIF_MIN_UNIT_SHARE: f64 = 0.1;
/// A further unit must differ from every unit already reported, of the same length, in at
/// least this fraction of its bases (over all rotations). CAG/CAA (a third) and AAAAG/AAGGG
/// (two fifths) are distinct units; the copies of a VNTR unit that differ in a few bases,
/// as at MUC1 (60-mer, ~5 differences) or TYMS, are variants of one.
const MOTIF_MIN_UNIT_DIFFERENCE: f64 = 0.2;
/// A further stretch, found in the sequence left on either side of the ones already used
/// (another period, such as a CT run beside a CCTT repeat), must score at least this
/// fraction of the first stretch.
const MOTIF_MIN_STRETCH_SHARE: f64 = 0.2;

/// The repeat units of an allele sequence, dominant first and joined by `/`
/// (e.g. `AAGGG/AAAAG`), or None when it holds no tandem repeat of a unit up to
/// `MOTIF_MAX_PERIOD` (see `MOTIF_MIN_SCORE`).
///
/// Per period, the best-scoring stretch where the sequence matches itself shifted by that
/// period is found (see `MOTIF_MISMATCH_PENALTY`); the period is the smallest one whose
/// stretch scores close to the best (see `MOTIF_PERIOD_MARGIN`), or a divisor of it that
/// scores reasonably (see `MOTIF_DIVISOR_MARGIN`). Codon repeats with a degenerate base,
/// such as polyalanine's GCN, are reported as their most common codon.
///
/// Within that stretch, units are counted with all rotations together (CAG, AGC and GCA
/// are one repeat read in different phases). A motif switch at the same period, as in
/// (AAAAG)n(AAGGG)m, stays in one stretch, since both halves match themselves shifted by
/// five, so every unit covering enough of it is reported (see `MOTIF_MIN_UNIT_SHARE`). The
/// sequence on either side of the stretch is then searched the same way for a repeat of
/// another period (see `MOTIF_MIN_STRETCH_SHARE`). Units are listed from the strongest
/// stretch to the weakest, and within a stretch from the most to the least frequent, at
/// most `MOTIF_MAX_UNITS`; a variant copy of a unit already listed is not a unit of its own
/// (see `MOTIF_MIN_UNIT_DIFFERENCE`).
///
/// Each unit is written in the rotation most common in the reference's frame, i.e. in
/// period-sized chunks from the start of the reference's own stretch for that period, when
/// the reference holds the same unit, else in the allele's. Anchoring on the reference makes
/// every allele and sample at a locus report the same rotation: HTT gives the same unit for
/// both alleles, rather than CAG for one and AGC for the other.
pub fn allele_motifs(seq: &str, reference: &str) -> Option<String> {
    let seq = seq.as_bytes().to_ascii_uppercase();
    let reference = reference.as_bytes().to_ascii_uppercase();
    // (windows covered, unit as its least rotation, period, the stretch it was found in), in
    // the order found: the strongest stretch first, most frequent unit first within one
    let mut units: Vec<(usize, Vec<u8>, usize, Range<usize>)> = Vec::new();
    let mut segments: Vec<Range<usize>> = Vec::new();
    segments.push(0..seq.len());
    let mut first_score = None;
    while units.len() < MOTIF_MAX_UNITS {
        // the best repeat among the segments not yet used
        let Some((index, score, period, stretch)) = segments
            .iter()
            .enumerate()
            .filter_map(|(i, segment)| {
                let (score, period, stretch) = dominant_period(&seq[segment.clone()])?;
                let shifted = segment.start + stretch.start..segment.start + stretch.end;
                Some((i, score, period, shifted))
            })
            .max_by_key(|(i, score, ..)| (*score, std::cmp::Reverse(*i)))
        else {
            break;
        };
        match first_score {
            None => first_score = Some(score),
            Some(first) if (score as f64) < first as f64 * MOTIF_MIN_STRETCH_SHARE => break,
            Some(_) => {}
        }
        let windows = stretch.len() + 1 - period;
        let classes = unit_counts(&seq[stretch.clone()], period);
        let dominant = classes.iter().map(|(_, count)| *count).max().unwrap_or(0);
        for (unit, count) in classes {
            // a unit that is itself a repeat, such as CTCT for a CT run that also matches
            // itself at the stretch's period of four, is its primitive unit
            let unit = primitive_unit(unit);
            let period = unit.len();
            let enough = count == dominant
                || (count as f64 >= windows as f64 * MOTIF_MIN_UNIT_SHARE && count >= 2 * period);
            if !enough
                || units
                    .iter()
                    .any(|(_, known, ..)| copy_variant(known, &unit))
            {
                continue;
            }
            units.push((count, unit, period, stretch.clone()));
        }
        let segment = segments.swap_remove(index);
        segments.extend(
            [segment.start..stretch.start, stretch.end..segment.end]
                .into_iter()
                .filter(|part| part.len() >= 2),
        );
    }
    if units.is_empty() {
        return None;
    }
    units.truncate(MOTIF_MAX_UNITS);
    let written: Vec<String> = units
        .into_iter()
        .map(|(_, unit, period, stretch)| {
            let reference_stretch = match reference.len() > period {
                true => best_stretch(&reference, period).1,
                false => 0..0,
            };
            let rotation = rotation_in_frame(&reference[reference_stretch], period, &unit)
                .or_else(|| rotation_in_frame(&seq[stretch], period, &unit))
                .unwrap_or(unit);
            String::from_utf8_lossy(&rotation).into_owned()
        })
        .collect();
    Some(written.join("/"))
}

/// (score, period, stretch) of the strongest tandem repeat in `seq`, see `allele_motifs`.
fn dominant_period(seq: &[u8]) -> Option<(i32, usize, Range<usize>)> {
    let max_period = MOTIF_MAX_PERIOD.min(seq.len() / 2);
    let stretches: Vec<(i32, Range<usize>)> =
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
        .filter(|d| candidate.is_multiple_of(*d))
        .find(|d| {
            let score = stretches[d - 1].0;
            score >= MOTIF_MIN_SCORE && score as f64 >= candidate_score * MOTIF_DIVISOR_MARGIN
        })
        .unwrap_or(candidate);
    let (score, stretch) = stretches[period - 1].clone();
    Some((score, period, stretch))
}

/// Whether `unit` is a variant copy of `known` rather than a unit of its own, see
/// `MOTIF_MIN_UNIT_DIFFERENCE`. Units of different lengths are always distinct.
fn copy_variant(known: &[u8], unit: &[u8]) -> bool {
    if known.len() != unit.len() {
        return false;
    }
    let n = unit.len();
    let closest = (0..n)
        .map(|r| (0..n).filter(|&i| known[i] != unit[(i + r) % n]).count())
        .min()
        .unwrap_or(0);
    (closest as f64) < MOTIF_MIN_UNIT_DIFFERENCE * n as f64
}

/// The shortest unit that repeats to make `unit` (CT for CTCT), as its least rotation.
fn primitive_unit(unit: Vec<u8>) -> Vec<u8> {
    let n = unit.len();
    match (1..n).find(|&d| n.is_multiple_of(d) && (d..n).all(|i| unit[i] == unit[i - d])) {
        Some(d) => least_rotation(&unit[..d]),
        None => unit,
    }
}

/// Windows of `seq` per unit, all rotations counted together, most frequent first (the
/// unit itself breaking ties, so the order is deterministic).
fn unit_counts(seq: &[u8], period: usize) -> Vec<(Vec<u8>, usize)> {
    // count raw k-mers first: far fewer distinct ones than windows in a repeat, so each
    // distinct k-mer is rotated once rather than once per position
    let mut raw: HashMap<&[u8], usize> = HashMap::new();
    for window in seq.windows(period) {
        *raw.entry(window).or_default() += 1;
    }
    let mut classes: BTreeMap<Vec<u8>, usize> = BTreeMap::new();
    for (kmer, count) in raw {
        *classes.entry(least_rotation(kmer)).or_default() += count;
    }
    let mut classes: Vec<(Vec<u8>, usize)> = classes.into_iter().collect();
    classes.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
    classes
}

/// The best-scoring stretch where `seq` matches itself shifted by `period` (maximum
/// subarray, Kadane), as its score and the bases it covers, both copies included.
fn best_stretch(seq: &[u8], period: usize) -> (i32, Range<usize>) {
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
    fn test_allele_motifs_period_and_phase() {
        let m = |seq: &str| allele_motifs(seq, "");
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
        // a motif switch at the same period reports both units, the larger first
        let mixed = format!("{}{}", "AAAAG".repeat(10), "AAGGG".repeat(40));
        assert_eq!(m(&mixed).as_deref(), Some("AAGGG/AAAAG"));
        // a unit covering a sliver of the allele is not reported
        let sliver = format!("{}{}", "AAAAG".repeat(2), "AAGGG".repeat(40));
        assert_eq!(m(&sliver).as_deref(), Some("AAGGG"));
        assert_eq!(m("cagcagcagcagcag").as_deref(), Some("CAG"));
        // interruptions are reported as a unit of their own, not as the multiple that spans
        // them: CAG with a CAA every sixth codon, or every third
        assert_eq!(m(&format!("{}CAA", "CAG".repeat(5)).repeat(3)).as_deref(), Some("CAG/CAA"));
        assert_eq!(m(&"CAGCAGCAA".repeat(4)).as_deref(), Some("CAG/CAA"));
        assert_eq!(m(&"GGGGCC".repeat(8)).as_deref(), Some("GGGGCC"));
        // the repeat inside flanking sequence, as in a padded catalog interval
        let flanked = format!("TACAGAAGTTTTTCTAGTATTGC{}ATAAGTTAAATAAAATTAAAT", "CAAA".repeat(8));
        assert_eq!(m(&flanked).as_deref(), Some("CAAA"));
    }

    #[test]
    fn test_allele_motifs_takes_the_phase_of_the_reference() {
        let reference = "CAG".repeat(19);
        // an allele that starts out of phase still reports the reference's rotation
        assert_eq!(
            allele_motifs(&format!("AG{}", "CAG".repeat(30)), &reference).as_deref(),
            Some("CAG")
        );
        // RFC1: an expanded AAAAG allele read in another phase is still AAAAG
        let rfc1 = "AAAAG".repeat(11);
        assert_eq!(
            allele_motifs(&format!("AAG{}", "AAAAG".repeat(100)), &rfc1).as_deref(),
            Some("AAAAG")
        );
        // a unit absent from the reference keeps the allele's own frame
        assert_eq!(allele_motifs(&"AAGGG".repeat(100), &rfc1).as_deref(), Some("AAGGG"));
    }

    #[test]
    fn test_least_rotation() {
        assert_eq!(least_rotation(b"CAG"), b"AGC");
        assert_eq!(least_rotation(b"GAAAA"), b"AAAAG");
        assert_eq!(least_rotation(b"ABAB"), b"ABAB");
        assert_eq!(least_rotation(b"A"), b"A");
        assert_eq!(primitive_unit(b"CTCT".to_vec()), b"CT");
        assert_eq!(primitive_unit(b"TCTC".to_vec()), b"CT");
        assert_eq!(primitive_unit(b"CCTT".to_vec()), b"CCTT");
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
    fn test_allele_motifs_reports_a_switch_to_another_period() {
        // a CT run beside a CCTT repeat, as at GOLGA8A; three distinct units at most
        let golga8a = format!("{}{}", "CT".repeat(15), "CCTT".repeat(30));
        let units = |seq: &str| -> Vec<Vec<u8>> {
            allele_motifs(seq, "")
                .unwrap()
                .split('/')
                .map(|u| least_rotation(u.as_bytes()))
                .collect()
        };
        assert_eq!(units(&golga8a), [b"CCTT".to_vec(), b"CT".to_vec()]);
        let four = format!(
            "{}{}{}{}",
            "CT".repeat(15),
            "CCTT".repeat(30),
            "GGA".repeat(20),
            "CA".repeat(12)
        );
        assert_eq!(allele_motifs(&four, "").map(|m| m.split('/').count()), Some(3));
        // a VNTR whose copies differ in a few bases is one unit
        let unit = "CCCGGGGCCGGCCTGGTGTCCGGGGCCGAGGTGACACCGTGGGCTGGGGGGGCGGTGGAG";
        let variant = unit
            .replacen("GGTGTC", "GCTCTC", 1)
            .replacen("GGGGGGG", "GCGGGCG", 1);
        let vntr = format!("{}{}{}{}", unit, variant, unit, variant);
        assert_eq!(allele_motifs(&vntr, "").map(|m| m.split('/').count()), Some(1));
        // errors in a pure repeat do not add units
        let mut noisy = "CAG".repeat(60).into_bytes();
        for i in (7..noisy.len()).step_by(23) {
            noisy[i] = b'T';
        }
        assert_eq!(allele_motifs(std::str::from_utf8(&noisy).unwrap(), "").as_deref(), Some("CAG"));
    }

    #[test]
    fn test_allele_motifs_none_when_not_a_repeat() {
        assert_eq!(allele_motifs("", ""), None);
        assert_eq!(allele_motifs(".", "."), None);
        assert_eq!(allele_motifs("ACGTTGCAAGCTTAGCCATGAGTCAGTACGATCGAATCGGCTA", ""), None);
    }

    #[test]
    #[ignore]
    fn test_create_motif() {
        assert_eq!(create_motif("CAGCAGCAGCAGCGGCGGCGGCAGCAGCAG"), "(CAG)4(CGG)3(CAG)3");
    }
}
