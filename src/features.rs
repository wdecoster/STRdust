//! Sequence feature vectors for the optional DBSCAN phasing strategy.
//!
//! Each insertion is mapped to a point `[length_weight * normalized_length,
//! canonical_kmer_freq_0, ...]`. The k-mer frequencies are length-invariant:
//! two reads of the same motif have near-identical k-mer vectors regardless of
//! how long the repeat is, which is what lets a length-variable expansion
//! cluster together (unlike the length-dominated Levenshtein distance used by
//! the default Ward strategy). Lifted from the `trout` cohort-outlier tool.
//!
//! Canonical k-mers fold cyclic rotations together (CAG ≡ AGC ≡ GCA) so the
//! representation is agnostic to where in the motif a read happens to start.

use std::collections::BTreeMap;

/// Bounds for `-k auto` period detection (k=1 is too crude; 4^6 already gives
/// ~700 canonical features, beyond which DBSCAN distances become noise).
pub const AUTO_K_MIN: usize = 2;
pub const AUTO_K_MAX: usize = 6;
/// Fallback when no clean period is detectable: trinucleotide composition has
/// the most signal among the common pathogenic motif lengths.
pub const AUTO_K_DEFAULT: usize = 3;
/// Match-rate threshold for `detect_period` (~35% mismatch tolerance).
const PERIOD_MATCH_THRESHOLD: f64 = 0.65;

/// Detect the dominant tandem-repeat period of `seq` by counting positions
/// where `seq[i] == seq[i+p]`. Returns the smallest p in `AUTO_K_MIN..=k_max`
/// clearing `PERIOD_MATCH_THRESHOLD` (the fundamental period), or None.
pub fn detect_period(seq: &[u8], k_max: usize) -> Option<usize> {
    if seq.len() < 2 * k_max {
        return None;
    }
    for p in AUTO_K_MIN..=k_max {
        let total = seq.len() - p;
        let matches = (0..total).filter(|&i| seq[i] == seq[i + p]).count();
        if (matches as f64) / (total as f64) >= PERIOD_MATCH_THRESHOLD {
            return Some(p);
        }
    }
    None
}

/// Precomputed mapping from every raw k-mer index to its canonical (lex-first
/// cyclic rotation) compact index.
pub struct KmerTable {
    /// raw_idx → compact canonical index (length: 4^k)
    pub compact: Vec<usize>,
    /// number of canonical classes
    pub n: usize,
}

impl KmerTable {
    pub fn new(k: usize) -> Self {
        let num_raw = 4usize.pow(k as u32);

        let decode = |idx: usize| -> Vec<usize> {
            (0..k)
                .map(|pos| (idx / 4usize.pow((k - 1 - pos) as u32)) % 4)
                .collect()
        };
        let encode = |bases: &[usize]| -> usize { bases.iter().fold(0, |acc, &b| acc * 4 + b) };

        // For each raw kmer, the raw index of its canonical (min) rotation.
        let canonical_raw: Vec<usize> = (0..num_raw)
            .map(|idx| {
                let bases = decode(idx);
                (0..k)
                    .map(|r| {
                        let rotated: Vec<usize> = (0..k).map(|i| bases[(i + r) % k]).collect();
                        encode(&rotated)
                    })
                    .min()
                    .unwrap_or(idx)
            })
            .collect();

        // Assign compact indices to canonical raws in sorted order (deterministic).
        let mut seen: BTreeMap<usize, usize> = BTreeMap::new();
        let mut next = 0usize;
        let compact: Vec<usize> = canonical_raw
            .iter()
            .map(|&cr| {
                *seen.entry(cr).or_insert_with(|| {
                    let i = next;
                    next += 1;
                    i
                })
            })
            .collect();

        KmerTable { compact, n: next }
    }
}

/// Build a feature vector: `[length_weight * normalized_length, canonical_kmer_freqs...]`.
pub fn build_feature_vector(
    seq: &str,
    table: &KmerTable,
    k: usize,
    max_len: usize,
    length_weight: f64,
) -> Vec<f64> {
    let normalized_len = if max_len > 0 {
        seq.len() as f64 / max_len as f64
    } else {
        0.0
    };
    let mut features = vec![normalized_len * length_weight];
    features.extend(kmer_frequencies(seq, table, k));
    features
}

fn kmer_frequencies(seq: &str, table: &KmerTable, k: usize) -> Vec<f64> {
    let seq_bytes = seq.as_bytes();
    if seq_bytes.len() < k {
        return vec![0.0; table.n];
    }
    let mut counts = vec![0u32; table.n];
    let total = seq_bytes.len() - k + 1;
    for window in seq_bytes.windows(k) {
        if let Some(raw) = kmer_raw_index(window) {
            counts[table.compact[raw]] += 1;
        }
    }
    let total_f = total as f64;
    counts.iter().map(|&c| c as f64 / total_f).collect()
}

fn kmer_raw_index(kmer: &[u8]) -> Option<usize> {
    let mut idx = 0usize;
    for &b in kmer {
        idx *= 4;
        idx += match b {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => return None,
        };
    }
    Some(idx)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kmer_frequencies_sum_to_one() {
        let table = KmerTable::new(2);
        let v = build_feature_vector("ACGTACGT", &table, 2, 8, 1.0);
        let sum: f64 = v[1..].iter().sum();
        assert!((sum - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_canonical_rotations_merged() {
        let table = KmerTable::new(3);
        // CAG, AGC, GCA all canonicalize to the same class
        let cag = kmer_raw_index(b"CAG").unwrap();
        let agc = kmer_raw_index(b"AGC").unwrap();
        let gca = kmer_raw_index(b"GCA").unwrap();
        assert_eq!(table.compact[cag], table.compact[agc]);
        assert_eq!(table.compact[cag], table.compact[gca]);
    }

    #[test]
    fn test_detect_period_cag() {
        assert_eq!(detect_period(b"CAGCAGCAGCAGCAGCAG", 6), Some(3));
    }

    #[test]
    fn test_length_invariance_of_composition() {
        // Same motif, very different lengths -> near-identical kmer part.
        let table = KmerTable::new(3);
        let short = "CAG".repeat(10);
        let long = "CAG".repeat(200);
        let vs = build_feature_vector(&short, &table, 3, long.len(), 1.0);
        let vl = build_feature_vector(&long, &table, 3, long.len(), 1.0);
        let comp_dist: f64 = vs[1..]
            .iter()
            .zip(vl[1..].iter())
            .map(|(a, b)| (a - b).powi(2))
            .sum();
        assert!(comp_dist < 1e-3, "composition should be length-invariant");
    }
}
