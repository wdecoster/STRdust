use crate::utils::levenshtein;
use kodama::{Method, linkage};
use log::{Level, debug, error, log_enabled};
use std::{cmp::max, collections::HashMap};

pub struct SplitSequences {
    pub hap1: Vec<String>,
    pub hap2: Option<Vec<String>>,
    pub flag: Option<String>,
    pub outliers: Option<Vec<String>>,
    /// Number of clusters found, only set by the DBSCAN strategy (None for Ward).
    /// Surfaced in the VCF as NCLUSTERS so complex multi-population loci can be flagged.
    pub n_clusters: Option<usize>,
}

pub fn split(
    insertions: &Vec<String>,
    repeat: &crate::repeats::RepeatInterval,
    check_outliers: bool,
    min_haplotype_fraction: f32,
    support: usize,
) -> SplitSequences {
    // the insertions are from an unphased experiment
    // and should be split in one (if homozygous) or two haplotypes
    // this is based on the length of the insertion
    // as well as the sequence composition
    // the sequence composition has to be agnostic to the motif length
    // and will use the levenshtein distance
    // this is inspired by the TRGT paper

    // Create a condensed (upper triangle) distance matrix
    // Create a cache for sequence lengths to avoid repeated calls to .len()
    let seq_lengths: Vec<usize> = insertions.iter().map(|s| s.len()).collect();

    // Early optimization: sort sequences by length to potentially exit Levenshtein calculation early
    let mut condensed = Vec::with_capacity((insertions.len() * (insertions.len() - 1)) / 2);

    // Using rayon for parallel processing
    use rayon::prelude::*;

    // Generate all pairs (row, col) that need processing
    let pairs: Vec<(usize, usize)> = (0..insertions.len() - 1)
        .flat_map(|row| ((row + 1)..insertions.len()).map(move |col| (row, col)))
        .collect();

    // First calculate the median length of all sequences
    let median_length = {
        let mut lengths = seq_lengths.clone();
        lengths.sort_unstable();
        if lengths.is_empty() {
            0 // Handle empty case, though this shouldn't happen
        } else if lengths.len().is_multiple_of(2) {
            (lengths[lengths.len() / 2] + lengths[lengths.len() / 2 - 1]) / 2
        } else {
            lengths[lengths.len() / 2]
        }
    };

    // Calculate distances in parallel with length-based penalties
    let distances: Vec<f32> = pairs
        .par_iter()
        .map(|(row, col)| {
            let len1 = seq_lengths[*row];
            let len2 = seq_lengths[*col];

            // Calculate basic length difference
            let len_diff = (len1 as i32 - len2 as i32).unsigned_abs() as usize;

            // For very different sequences, use length difference but with moderate scaling
            if len_diff > 100 {
                return len_diff as f32 * 1.2; // Apply a modest 20% increase instead of more aggressive scaling
            }

            // For more similar sequences, calculate Levenshtein distance
            let base_distance = levenshtein(&insertions[*row], &insertions[*col]) as f32;

            // Apply a gentle penalty based on relative length difference compared to median
            let longer_len = len1.max(len2);

            // Only apply additional penalty when:
            // 1. There's a meaningful but not extreme length difference
            // 2. One sequence is notably larger than median
            if len_diff > 20
                && len_diff <= 100
                && (longer_len as f32) > (median_length as f32 * 1.3)
            {
                // Calculate how much the longer sequence exceeds median (as a ratio)
                let expansion_factor = (longer_len as f32 / median_length as f32) - 1.0;

                // Apply a moderate, logarithmic penalty that increases more slowly
                // This avoids overly aggressive penalties while still separating different haplotypes
                let penalty_factor = 1.0 + (expansion_factor * 0.5).min(1.0);

                return base_distance * penalty_factor;
            }

            // For sequences with similar lengths or not significantly expanded, just use Levenshtein
            base_distance
        })
        .collect();

    condensed.extend(distances);
    let dend = linkage(&mut condensed, insertions.len(), Method::Ward);

    // maybe these hashmaps could be replaced by some tree-like structure
    // create a hashmap to store labels of the clusters and their subclusters
    let mut cluster_to_subclusters = HashMap::new();
    // create a hashmap to store the size of the clusters
    let mut clusters_to_size = HashMap::new();
    // create a hashmap to store the dissimilarity of the clusters
    let mut clusters_to_dissimilarity = HashMap::new();
    // create a hashmaps to store the subclusters and their parent cluster
    let mut subcluster_to_cluster = HashMap::new();
    // create a vector to store the clusters
    let mut clusters = vec![];
    // create a vector with candidate haplotype-clusters
    let mut haplotype_clusters = vec![];
    // clusters have to represent at least min_haplotype_fraction of the reads, but never less than 1
    // the fraction-based threshold is capped by `support` so that copy-number gains (which inflate
    // the total read count, and thus the fraction-based threshold) cannot swallow a genuine minority
    // allele that is supported by at least `support` reads
    let min_cluster_size =
        max(((insertions.len() as f32 * min_haplotype_fraction) as usize).min(support), 1);
    debug!(
        "{repeat}: Minimum cluster size: {} reads (min of {}% of {} total and support {})",
        min_cluster_size,
        min_haplotype_fraction * 100.0,
        insertions.len(),
        support
    );

    for (index, step) in dend.steps().iter().enumerate() {
        // insert the new label with the clusters it contains
        // at every step of the dendrogram the label of the cluster is N + i
        // where N is the number of observations and i is the index of the step
        // the cluster hashmap is updated with the new label and the clusters it contains
        let label = index + insertions.len();
        cluster_to_subclusters.insert(label, (step.cluster1, step.cluster2));
        clusters_to_size.insert(label, step.size);
        clusters_to_dissimilarity.insert(label, step.dissimilarity);
        subcluster_to_cluster.insert(step.cluster1, label);
        subcluster_to_cluster.insert(step.cluster2, label);
        // clusters are added to a vector, as a tuple with their size
        clusters.push((index + insertions.len(), step.size));
    }
    // sort the clusters by size
    clusters.sort_by_key(|x| x.1);
    clusters.reverse();

    // when debugging, the code below can be used to print the tree
    if log_enabled!(Level::Debug) {
        for (cluster, _size) in clusters.iter() {
            let subclusters = cluster_to_subclusters.get(cluster).unwrap();

            let seq1 = if subclusters.0 < insertions.len() {
                insertions[subclusters.0].clone()
            } else {
                "".to_string()
            };
            let seq2 = if subclusters.1 < insertions.len() {
                insertions[subclusters.1].clone()
            } else {
                "".to_string()
            };

            debug!(
                "{repeat}: Node {cluster} with dissimilarity {} and children {} [{seq1}] and {} [{seq2}]",
                clusters_to_dissimilarity.get(cluster).unwrap(),
                subclusters.0,
                subclusters.1
            );
        }
    }
    // figure out which cluster should be considered root(s)
    let roots = find_roots(
        clusters[0].0,
        &cluster_to_subclusters,
        &clusters_to_size,
        &clusters_to_dissimilarity,
        &min_cluster_size,
    );
    debug!("{repeat}: Roots for this tree: {:?}", roots);

    // if a parent cluster has been seen we will ignore all children thereof
    let mut large_cluster_seen = vec![];
    for (cluster, size) in clusters.iter() {
        // have to ignore the biggest cluster - the one with all reads, as well as clusters that are too small
        // in the case of an outlier the next cluster will be the largest and has to be considered the root
        if !roots.contains(cluster) && size > &min_cluster_size {
            // for every cluster, find the parent cluster
            let parent = if !subcluster_to_cluster.contains_key(cluster) {
                // if the cluster is not in the hashmap, it is a root
                // and therefore its parent is the cluster itself
                cluster
            } else {
                // otherwise the parent is the cluster in the hashmap
                subcluster_to_cluster.get(cluster).unwrap()
            };
            // if the parent cluster has already been seen for this large cluster, we ignore it
            // as such we only get sufficiently large independent clusters
            if !large_cluster_seen.contains(parent) {
                haplotype_clusters.push(*cluster);
                debug!("{repeat}: Adding cluster {} to candidate haplotype clusters", cluster);
            }
            large_cluster_seen.push(*cluster);
        }
    }
    match haplotype_clusters.len() {
        0 => {
            debug!(
                "{repeat}: No haplotype clusters found! Treating this as homozygous, but here could be dragons"
            );
            SplitSequences {
                hap1: insertions.clone(),
                hap2: None,
                flag: Some("CLUSTERFAILURE".to_string()),
                // looking for outlier lengths that could be a poorly covered allele
                outliers: if check_outliers {
                    find_outliers(insertions, None)
                } else {
                    None
                },
                n_clusters: None,
            }
        }
        1 => {
            // if there is only one haplotype cluster, and the locus is considered homozygous
            // all insertions are returned as the first haplotype with None as the second haplotype
            // this ignores roots and noise insertions, but identifying those is rather problematic in the homozygous case
            // and the end result is often that we lose too many reads as false-positive roots
            // I assume the poa consensus will deal with noise insertions
            debug!("{repeat}: Only one haplotype cluster found");
            SplitSequences {
                hap1: insertions.clone(),
                hap2: None,
                flag: None,
                outliers: if check_outliers {
                    find_outliers(insertions, None)
                } else {
                    None
                },
                n_clusters: None,
            }
        }
        2 => {
            debug!("{repeat}: Found two haplotype clusters");
            let hap1_refs =
                find_cluster_members(&haplotype_clusters[0], &cluster_to_subclusters, insertions);
            let hap1: Vec<String> = hap1_refs.iter().map(|s| (*s).clone()).collect();
            let hap2_refs =
                find_cluster_members(&haplotype_clusters[1], &cluster_to_subclusters, insertions);
            let hap2: Vec<String> = hap2_refs.iter().map(|s| (*s).clone()).collect();
            let larger_median = if check_outliers {
                let hap1_median = find_median(&hap1);
                let hap2_median = find_median(&hap2);
                Some(max(hap1_median, hap2_median))
            } else {
                None
            };
            SplitSequences {
                hap1,
                hap2: Some(hap2),
                flag: None,
                outliers: if check_outliers {
                    find_outliers(insertions, larger_median)
                } else {
                    None
                },
                n_clusters: None,
            }
        }
        _ => {
            error!("{repeat}: Found more than two haplotype clusters. This shouldn't happen.");
            panic!();
        }
    }
}

/// Alternative phasing strategy to [`split`] (Ward over length-weighted
/// Levenshtein). Builds k-mer composition feature vectors (length-invariant)
/// and clusters them with DBSCAN. Noise points and any clusters beyond the two
/// largest become outliers. Returns `n_clusters` so complex multi-population
/// loci can be flagged in the VCF.
pub fn split_dbscan(
    insertions: &[String],
    repeat: &crate::repeats::RepeatInterval,
    check_outliers: bool,
    support: usize,
    eps: f64,
    length_weight: f64,
) -> SplitSequences {
    use crate::features::{self, KmerTable};

    // Auto-detect the repeat period (k) from the median-length insertion.
    let mut by_len: Vec<&String> = insertions.iter().collect();
    by_len.sort_by_key(|s| s.len());
    let median_seq = by_len[by_len.len() / 2];
    let k = features::detect_period(median_seq.as_bytes(), features::AUTO_K_MAX)
        .unwrap_or(features::AUTO_K_DEFAULT);
    let table = KmerTable::new(k);
    let max_len = insertions.iter().map(|s| s.len()).max().unwrap_or(0);

    let points: Vec<Vec<f64>> = insertions
        .iter()
        .map(|s| features::build_feature_vector(s, &table, k, max_len, length_weight))
        .collect();

    // min_samples = support, so every DBSCAN cluster is backed by >= support reads.
    let labels = crate::dbscan::cluster(&points, eps, support);

    let mut clusters: HashMap<usize, Vec<String>> = HashMap::new();
    let mut noise: Vec<String> = Vec::new();
    for (seq, label) in insertions.iter().zip(labels.iter()) {
        match label {
            Some(cid) => clusters.entry(*cid).or_default().push(seq.clone()),
            None => noise.push(seq.clone()),
        }
    }
    let n_clusters = Some(clusters.len());
    debug!(
        "{repeat}: DBSCAN (k={k}, eps={eps}, min_samples={support}) found {} clusters and {} noise reads",
        clusters.len(),
        noise.len()
    );

    // Sort clusters by size (desc); tie-break on total bp then id for determinism.
    let mut sized: Vec<(usize, Vec<String>)> = clusters.into_iter().collect();
    sized.sort_by(|a, b| {
        b.1.len()
            .cmp(&a.1.len())
            .then_with(|| {
                b.1.iter()
                    .map(String::len)
                    .sum::<usize>()
                    .cmp(&a.1.iter().map(String::len).sum::<usize>())
            })
            .then_with(|| a.0.cmp(&b.0))
    });

    let outliers_opt = |extra: Vec<String>| -> Option<Vec<String>> {
        check_outliers.then_some(extra).filter(|e| !e.is_empty())
    };

    match sized.len() {
        0 => {
            debug!("{repeat}: DBSCAN found no clusters; treating as homozygous (CLUSTERFAILURE)");
            SplitSequences {
                hap1: insertions.to_vec(),
                hap2: None,
                flag: Some("CLUSTERFAILURE".to_string()),
                outliers: None,
                n_clusters,
            }
        }
        1 => {
            let hap1 = sized.into_iter().next().unwrap().1;
            SplitSequences {
                hap1,
                hap2: None,
                flag: None,
                outliers: outliers_opt(noise),
                n_clusters,
            }
        }
        _ => {
            let mut it = sized.into_iter();
            let hap1 = it.next().unwrap().1;
            let hap2 = it.next().unwrap().1;
            // any clusters beyond the two largest join the noise as outliers
            let mut extra = noise;
            for (_, c) in it {
                extra.extend(c);
            }
            SplitSequences {
                hap1,
                hap2: Some(hap2),
                flag: None,
                outliers: outliers_opt(extra),
                n_clusters,
            }
        }
    }
}

fn find_roots(
    top_root: usize,
    cluster_to_subclusters: &HashMap<usize, (usize, usize)>,
    clusters_to_size: &HashMap<usize, usize>,
    clusters_to_dissimilarity: &HashMap<usize, f32>,
    min_cluster_size: &usize,
) -> Vec<usize> {
    // find nodes that qualify as roots
    // this includes the top most node, for which the size is equal to the number of sequences
    // as well as its children nodes that are the sibling of a too small node
    let mut roots = vec![top_root];
    // if the cluster is not in the hashmap, then move on
    if !cluster_to_subclusters.contains_key(&top_root) {
        return roots;
    }
    let (child1, child2) = cluster_to_subclusters.get(&top_root).unwrap();
    let size1 = clusters_to_size.get(child1).unwrap_or(&0);
    let size2 = clusters_to_size.get(child2).unwrap_or(&0);
    if size1 > min_cluster_size && size2 > min_cluster_size {
        // if both clusters are sufficiently large we are done finding roots
        return roots;
    } else if clusters_to_dissimilarity.get(&top_root).unwrap() < &5.0 {
        // if one node is too small, but the difference is not large enough, there is no root that has to be ignored
        return vec![];
    } else {
        // if one of the clusters is too small and the dissimilarity is large, the other one is a root
        // in that case we have to recurse to find the roots of the other cluster
        if size1 > min_cluster_size {
            roots.extend(find_roots(
                *child1,
                cluster_to_subclusters,
                clusters_to_size,
                clusters_to_dissimilarity,
                min_cluster_size,
            ));
        } else {
            roots.extend(find_roots(
                *child2,
                cluster_to_subclusters,
                clusters_to_size,
                clusters_to_dissimilarity,
                min_cluster_size,
            ));
        }
    }
    roots
}

fn find_cluster_members<'a>(
    cluster: &usize,
    cluster_to_subclusters: &HashMap<usize, (usize, usize)>,
    insertions: &'a [String],
) -> Vec<&'a String> {
    let mut result = Vec::with_capacity(insertions.len() / 2);
    let mut to_process = vec![*cluster];

    while let Some(current) = to_process.pop() {
        if let Some(&(c1, c2)) = cluster_to_subclusters.get(&current) {
            if c1 < insertions.len() {
                result.push(&insertions[c1]);
            } else {
                to_process.push(c1);
            }

            if c2 < insertions.len() {
                result.push(&insertions[c2]);
            } else {
                to_process.push(c2);
            }
        }
    }

    result
}

fn find_median(seqs: &Vec<String>) -> usize {
    let mut lengths = vec![];
    for seq in seqs {
        lengths.push(seq.len());
    }
    // sort the lengths
    lengths.sort_unstable();
    // find the median length, but check if there is an even or odd number of lengths
    if lengths.len() % 2 == 0 {
        (lengths[lengths.len() / 2] + lengths[lengths.len() / 2 - 1]) / 2
    } else {
        lengths[lengths.len() / 2]
    }
}

fn find_outliers(seqs: &Vec<String>, larger_median: Option<usize>) -> Option<Vec<String>> {
    // this generates optional output only to be ran with --find-outliers
    // find outlier insertions that are much longer than the rest as these could be a poorly covered allele
    // this is based on the length of the insertion and an insertion is considered an outlier if it is more than twice the median length
    // in the case of a heterozygous allele the median is already calculated
    // and we will use the median of the larger allele to identify outliers
    let median_length = match larger_median {
        Some(median) => median,
        None => find_median(seqs),
    };

    let mut outliers = vec![];
    for seq in seqs.iter() {
        if seq.len() > median_length * 2 {
            outliers.push(seq.clone());
        }
    }
    debug!("Found {} outliers.", outliers.len());
    if outliers.is_empty() {
        None
    } else {
        Some(outliers)
    }
}

mod tests {
    #[allow(unused_imports)]
    use super::*;

    /// The 44 insertions of chr15:34419425-34419451 (GOLGA8A) from rr_NA99_170.log:
    /// 5 length-variable, CT/CCTT-rich expansion reads + 39 short TTTC reference reads.
    #[allow(dead_code)]
    fn golga8a_insertions() -> Vec<String> {
        let r43 = "TCTTTCTTTCTTTCCTTTCCTTTCCTTTCCTTTCCTTTCCTTCCTTCCCTTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTTTTTCTTTCTTTCTTT".to_string();
        let r2 = "TCTTTCTTTCTTCCTTTCCTTTCCTTTCCTTTCCTTTCCTTCCTTCCCTTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTTTCTTTCTTTCTTT".to_string();
        let r33 = "CTTTCTTTCCTTTCCTTTCCTTTCCTTTCCTTTCCTTCCTTCCCCCTGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAAAAAAAAAAAAAAAAAAAAAAAATCTCTCTCTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTTTCTTTCTTTCTTTCTTTCTTTC".to_string();
        let r8 = "TCTTTCTTTCTTTCCTTTCCTTTCCTTTCCTTTCCTTTCCTTCCTTCCCTTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTTTCTTCTTTCTTT".to_string();
        let r23 = "TCTTTCTTTCCTTTCCTTTCCTTTCCTTTCCTTTCCTTCCTTCCCTTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTTCTCTCTCTCTCTCTCTCTCTTTCTTTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTTTCTTTCTTTCTTT".to_string();

        let mut insertions = vec![r43, r2, r33, r8, r23];
        // 39 short reference-allele reads (lengths matching the logged distribution)
        for _ in 0..20 {
            insertions.push("TTTCTTTCTTTCTTTCTTTCTTTCTTT".to_string()); // 27
        }
        for _ in 0..7 {
            insertions.push("TTTCTTTCTTTCTTTCTTTCTTTCTTTT".to_string()); // 28
        }
        for _ in 0..2 {
            insertions.push("TTTCTTTCTTTCTTTCTTTCTTTCTT".to_string()); // 26
        }
        for _ in 0..2 {
            insertions.push("TTTCTTTCTTTCTTTCTTTCTTTCTTTTT".to_string()); // 29
        }
        insertions.push("CTTTTCTTTCTTTCTTTCTTTC".to_string()); // 22
        insertions.push("TTTCTTTCTTTCTTTCTTTCTTT".to_string()); // 23
        insertions.push("TTTTCTTTCTTTCTTTCTTTCTTT".to_string()); // 24
        insertions.push("TTTCTTTCTTTCTTTCTTTCTTTCTTTCTTT".to_string()); // 31
        insertions.push("AGTTTCTTTCTTTCTTTCTTTCTTTTTTTCTTT".to_string()); // 33
        insertions.push("TTTCTTTCTTTCTTTCTTTTTTTCTTTTTTTCTTT".to_string()); // 35
        for _ in 0..2 {
            insertions.push("TTTCTTTCTTTCTTTCTTTCTTTCTTTTTTTTTT".to_string()); // 34
        }
        assert_eq!(insertions.len(), 44);
        insertions
    }

    #[test]
    fn test_split_dbscan_captures_expansion() {
        // The DBSCAN strategy on the GOLGA8A data: with length-invariant k-mer
        // composition features, the length-variable expansion is recovered as a
        // distinct haplotype (rather than being fragmented by length and lost to
        // the outlier bin, as happens with the length-dominated Ward distance).
        let insertions = golga8a_insertions();
        let repeat = crate::repeats::RepeatInterval {
            chrom: "chr15".to_string(),
            start: 34419425,
            end: 34419451,
            created: None,
        };
        // default CLI params (eps=0.4, length_weight=0.3), min_samples = support = 2
        let s = split_dbscan(&insertions, &repeat, true, 2, 0.4, 0.3);
        let hap2 = s.hap2.expect("expected a heterozygous call");
        // exactly the reference + expansion alleles -> two clusters
        assert_eq!(s.n_clusters, Some(2));
        // one haplotype is the short reference, the other captures the expansion
        let med1 = find_median(&s.hap1);
        let med2 = find_median(&hap2);
        assert!(
            med1.min(med2) < 50 && med1.max(med2) > 100,
            "expected one reference (~27bp) and one expansion (>100bp) haplotype, got {med1} and {med2}"
        );
    }

    #[test]
    fn test_split_expansion_minority_allele() {
        // Regression test for chr15:34419425-34419451 from rr_NA99_170.log
        // A real, length-variable expansion is present on only 5/44 reads (the rest
        // are the short reference allele). The expansion reads must be recovered as a
        // haplotype and must NOT all be discarded as length outliers.
        let mut insertions = golga8a_insertions();

        use rand::seq::SliceRandom;
        let mut rng = rand::rng();
        insertions.shuffle(&mut rng);

        let splitseqs = split(
            &insertions,
            &crate::repeats::RepeatInterval {
                chrom: "chr15".to_string(),
                start: 34419425,
                end: 34419451,
                created: None,
            },
            true,
            0.1,
            2,
        );

        let hap2 = splitseqs
            .hap2
            .expect("expected a heterozygous call with two haplotypes");
        // One of the two haplotypes should capture the expansion (median length >> reference)
        let expansion_captured = find_median(&splitseqs.hap1) > 100 || find_median(&hap2) > 100;
        let n_outliers = splitseqs.outliers.as_ref().map_or(0, |o| o.len());
        assert!(
            expansion_captured,
            "expansion not captured as a haplotype; hap1 median {}, hap2 median {}, {} outliers",
            find_median(&splitseqs.hap1),
            find_median(&hap2),
            n_outliers
        );
    }

    #[test]
    fn test_split_length() {
        // test that the split function identifies two haplotype, mainly based on insertion length
        // using a CAG expansion with some SNV noise
        let mut insertions = vec![
            "CTGCAGCAGCAGCTGCAGCAGCAGCAGCAGCAGCAGTAGCAGCAGCAGCAGCTGCAGCAG".to_string(),
            "CAGCGGCAGCAGCAGCAGCAGCAGCAGCAGCAGCGGCAGCAGCAGCAGCAGCAGCAGCAG".to_string(),
            "CAGCTGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCTGCAGCAGGAGCATCAG".to_string(),
            "CAGCAGCAGCAGCTGCAGCAGCAGCAGCAG".to_string(),
            "CAGCCGCAGCAGCAGCAGCAGCAGCAGCAG".to_string(),
            "CAGCAGCGGCAGCAGCAGCAGCAGCAGCAG".to_string(),
        ];
        // shuffle the insertions in a random order
        use rand::seq::SliceRandom;
        let mut rng = rand::rng();
        insertions.shuffle(&mut rng);
        let splitseqs = split(
            &insertions,
            &crate::repeats::RepeatInterval {
                chrom: "chr7".to_string(),
                start: 154654404,
                end: 154654432,
                created: None,
            },
            false,
            0.1,
            3,
        );
        assert!(splitseqs.hap1.len() == splitseqs.hap2.unwrap().len());
        // check that all sequences in hap1 are the same length
        assert!(
            splitseqs
                .hap1
                .iter()
                .all(|x| x.len() == splitseqs.hap1[0].len())
        );
    }

    #[test]
    fn test_split_composition() {
        // test that the split function identifies two haplotype, based on sequence composition
        // using a AAGGG and AAAAG expansion with some SNV and length noise
        // this is the composition in the RFC1 expansion in CANVAS, with the former being pathogenic
        let mut insertions = vec![
            "AAGGGAAGGGAAGGGAAGGGAATGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGG".to_string(),
            "AAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGAAGGGAAGGG".to_string(),
            "AAGGGAAGGGAAGGGAAGGGAAGGGCAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGG".to_string(),
            "AAAAGAAAAGAAAAGAAAAGAAAAGAAAGGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAG".to_string(),
            "AAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAGGAAAAGAAAAGAAGAAAAGAAAAGAAAAG".to_string(),
            "AAAAGAAAAGAAAAGAAAAGAAAAAGAAAAGAAAAGAAAAGAAACGAAAAGAAAAGAAAAGAAAAG".to_string(),
        ];
        // the expected_haplotype are the first three insertions
        let mut expected_haplotype = vec![
            "AAGGGAAGGGAAGGGAAGGGAATGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGG".to_string(),
            "AAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGAAGGGAAGGG".to_string(),
            "AAGGGAAGGGAAGGGAAGGGAAGGGCAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGG".to_string(),
        ];
        // shuffle the insertions in a random order
        use rand::seq::SliceRandom;
        let mut rng = rand::rng();
        insertions.shuffle(&mut rng);
        let splitseqs = split(
            &insertions,
            &crate::repeats::RepeatInterval {
                chrom: "chr7".to_string(),
                start: 154654404,
                end: 154654432,
                created: None,
            },
            false,
            0.1,
            3,
        );
        let mut hap1 = splitseqs.hap1;
        let mut hap2 = splitseqs.hap2.unwrap();
        assert!(hap1.len() == hap2.len());
        // assert that either hap1 or hap2 is equal to the expected_haplotype
        // but first sort the haplotypes, since the order is not guaranteed
        hap1.sort();
        hap2.sort();
        expected_haplotype.sort();
        assert!(hap1 == expected_haplotype || hap2 == expected_haplotype);
    }

    #[test]
    fn test_split_with_homozygous() {
        // test how the split function behaves with one 'homozygous' haplotype,
        // using a CAG expansion with some SNV and length noise
        // this will still return two haplotypes, as there will always be a tree with two branches
        // if the sequences are truly similar enough then the consensus will be about the same
        // or at least within the edit distance cutoff
        let mut insertions = vec![
            "CTGCAGCAGCAGCTGCAGCAGCAGCAGCAGCAGCAGTAGCAGCAGCAGCAGCTGCAGCAG".to_string(),
            "CAGCGGCAGCAGCAGCAGCAGCGCAGCAGCAGCGGCAGCAGCAGCAGCAGCAGCAGCAG".to_string(),
            "CAGCTGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCTGCAGCAGGAGCATCAG".to_string(),
            "CGGCAGCAGCAGCTGCAGCAGCAGCAGCAGCAGCAGTAGCAGCAGCAGCAGCTGCAGCAG".to_string(),
            "CAGCGGCAGCAGCACCCAGCACCAGCAGCAGCAGCGGCAGCAGCAGCAGCAGCAGCAGCAG".to_string(),
            "CAGCTGCAGCAGCAGCAGCAGCAGCGGCCGCAGCAGCAGCAGCTGCAGCAGGAGCATCAG".to_string(),
        ];
        // shuffle the insertions in a random order
        use rand::seq::SliceRandom;
        let mut rng = rand::rng();
        insertions.shuffle(&mut rng);
        let splitseqs = split(
            &insertions,
            &crate::repeats::RepeatInterval {
                chrom: "chr7".to_string(),
                start: 154654404,
                end: 154654432,
                created: None,
            },
            false,
            0.1,
            3,
        );
        assert!(splitseqs.hap1.len() + splitseqs.hap2.unwrap().len() == insertions.len());
    }

    #[test]
    fn test_split_with_homozygous2() {
        // this is a real example of a homozygous one, where the tree is just one long trail of nodes
        // where every time a single new sequence joins the tree
        // this therefore created a problem in find_roots, as every node became a root
        let mut insertions = vec![
            "CAGCAGCAGCAGCA".to_string(),
            "CAACAGCAGCAGCAGCA".to_string(),
            "CAGCAGCAGCAGCAGCA".to_string(),
            "CAGCAGCAGCAGCAGCA".to_string(),
            "CAGCAGCAGCAGCAGCA".to_string(),
            "CAGCAGCAGCAGCAGCA".to_string(),
            "CAGCAGCAGCAGCAGCA".to_string(),
            "CAGCAGCAGCAGCAGCA".to_string(),
            "CAGCAGCAGCAGCAGCA".to_string(),
            "CAGCAGCAGCAGCAGCA".to_string(),
            "CAGCAGCAGCAGCAGCA".to_string(),
            "CAGCAGCAGCAGCAGCA".to_string(),
            "CAGCAGCAGCAGCAGCA".to_string(),
            "CAGCAGCAGCAGCAGCA".to_string(),
            "CAGCAGCAGCAGCAGCA".to_string(),
            "CAGCAGCAGCAGCAGCA".to_string(),
            "CAGCAGCAGCAGCAGCA".to_string(),
            "CAGCAGCAGCAGCAGCA".to_string(),
            "CAGCAGCAGCAGCAGCA".to_string(),
        ];
        // shuffle the insertions in a random order
        use rand::seq::SliceRandom;
        let mut rng = rand::rng();
        insertions.shuffle(&mut rng);
        let splitseqs = split(
            &insertions,
            &crate::repeats::RepeatInterval {
                chrom: "chr7".to_string(),
                start: 154654404,
                end: 154654432,
                created: None,
            },
            false,
            0.1,
            3,
        );
        assert!(splitseqs.hap2.is_none());
        println!("hap1: {:?}", splitseqs.hap1);
        assert!(splitseqs.hap1.len() == insertions.len());
    }

    #[test]
    fn test_split_with_noise() {
        // test that the split function identifies two haplotype, based on sequence composition
        // using a AAGGG and AAAAG expansion with some SNV and length noise
        // this is the composition in the RFC1 expansion in CANVAS, with the former being pathogenic
        // additionally, an outlier insertion is added that does not belong to any of the haplotypes
        // this outlier should be discarded, returning the same output as the test above
        let mut insertions = vec![
            "AAGGGAAGGGAAGGGAAGGGAATGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGG".to_string(),
            "AAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGAAGGGAAGGG".to_string(),
            "AAGGGAAGGGAAGGGAAGGGAAGGGCAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGG".to_string(),
            "AAAAGAAAAGAAAAGAAAAGAAAAGAAAGGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAG".to_string(),
            "AAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAGGAAAAGAAAAGAAGAAAAGAAAAGAAAAG".to_string(),
            "AAAAGAAAAGAAAAGAAAAGAAAAAGAAAAGAAAAGAAAAGAAACGAAAAGAAAAGAAAAGAAAAG".to_string(),
            "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT".to_string(),
        ];
        // the expected_haplotype are the first three insertions
        let mut expected_haplotype = vec![
            "AAGGGAAGGGAAGGGAAGGGAATGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGG".to_string(),
            "AAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGAAGGGAAGGG".to_string(),
            "AAGGGAAGGGAAGGGAAGGGAAGGGCAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGG".to_string(),
        ];
        // shuffle the insertions in a random order
        use rand::seq::SliceRandom;
        let mut rng = rand::rng();
        insertions.shuffle(&mut rng);
        let splitseqs = split(
            &insertions,
            &crate::repeats::RepeatInterval {
                chrom: "chr7".to_string(),
                start: 154654404,
                end: 154654432,
                created: None,
            },
            false,
            0.1,
            3,
        );
        let mut hap1 = splitseqs.hap1;
        let mut hap2 = splitseqs.hap2.unwrap();
        println!("{:?}", hap1);
        println!("{}", hap1.len());
        println!("{:?}", hap2);
        assert!(hap1.len() == hap2.len());
        // assert that either hap1 or hap2 is equal to the expected_haplotype
        // but first sort the haplotypes, since the order is not guaranteed
        hap1.sort();
        hap2.sort();
        expected_haplotype.sort();
        assert!(hap1 == expected_haplotype || hap2 == expected_haplotype);
    }

    #[test]
    fn test_split_with_two_times_noise() {
        // test that the split function identifies two haplotype, based on sequence composition
        // using a AAGGG and AAAAG expansion with some SNV and length noise
        // this is the composition in the RFC1 expansion in CANVAS, with the former being pathogenic
        // additionally, two outlier insertions are added that do not belong to any of the haplotypes
        // these outliers should be discarded
        // additional 'good' reads are added to increase the coverage
        let mut insertions = vec![
            "AAGGGAAGGGAAGGGAAGGGAATGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGG".to_string(),
            "AAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGAAGGGAAGGG".to_string(),
            "AAGGGAAGGGAAGGGATGGGAAGGGAAGGGAGGGGAAGGGAAGGAACAGGGGAAGGGAAGGGAAGGG".to_string(),
            "AAGGGAAGGGAAGGGAAGGGAAGGGCAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGG".to_string(),
            "AAGGGAAGGGAAGGGATGGGAAGGGCAGGGAAGGGAAGGGAAGGAAAAGGGAAGGGAAGGGAAGGG".to_string(),
            "AAGGGAAGAATGGGAAGGGAAGGGCAGGGAAGGGAAGGGGAAGGGAAGGGAAGGGAAGGGAAGGG".to_string(),
            "AAGGGAAGGGAAGGGATGGGAAGGGCAGGGAAGGGAAGGGAAGGAAAAGGGAAGGGAAGGGAAGGG".to_string(),
            "AAGGGAAGAATGGGAAGGGAAGGGCAGGGAAGGGAAGGGGAAGGGAAGGGAAGGGAAGGGAAGGG".to_string(),
            "AAGGGGGGAAGGGACGGGAAGGGCAGGGAAGGGAAGGGAAGGAAAAAGGGAAGGGAAGGGAAGGG".to_string(),
            "AAAAGAAAAGAAAAGAAAAGAAAAGAAAGGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAG".to_string(),
            "AAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAGGAAAAGAAAAGAAGAAAAGAAAAGAAAAG".to_string(),
            "AAAAGAAAAGAAAAGAAAAGAAAAAGAAAAGAAAAGAAAAGAAACGAAAAGAAAAGAAAAGAAAAG".to_string(),
            "AAAAAAGAAAGAAAGAAAAGAAAAAGGGAAAAGAAAAGAAAAGAAACGAAAAGAAAAGAAAAGAAAAAG".to_string(),
            "AAAAGGAAAGAAAGAAAAGAAAAAGGGAAAAGAAAAGAAAAGAAACGAAAAGAAAAGAAAAGAAAAAG".to_string(),
            "AAAAGAAGGGAAGAAAAGAAAAAGGGAAAAGAAAAGAAAAGAAACGAAAAGAAAAGAAAAGAAAAAG".to_string(),
            "AAAAAAGAAAGAAAGAAAAGAAAAAGGGAAAAGAAAAGAAAAGAAACGAAAAGAAAAGAAAAGAAAAAG".to_string(),
            "AAAAGGAAAGAAAGAAAAGAAAAAGGGAAAAGAAAAGAAAAGAAACGAAAAGAAAAGAAAAGAAAAAG".to_string(),
            "AAAAGAAGGGAAGAAAAGAAAAAGGGAAAAGAAAAGAAAAGAAACGAAAAGAAAAGAAAAGAAAAAG".to_string(),
            "GAAAGAAAGAAAAGAAAAAGAGGAAAAGGAAAGAAAAGAAACGAAAAGAAAAGAAAAGAAAAAG".to_string(),
            "AAAAGGAAAGGGAAGAAAAGAAAAAGAGGAAAAGAAAAGAAAAGAAACGGGAAAAGAAAAGAAGAAAAAG".to_string(),
            "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT".to_string(),
            "CCCCCCCCCCCCCCCCCCCCCC".to_string(),
        ];
        // the expected_haplotype are the first three insertions
        let mut expected_haplotype = vec![
            "AAGGGAAGGGAAGGGAAGGGAATGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGG".to_string(),
            "AAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGAAGGGAAGGG".to_string(),
            "AAGGGAAGGGAAGGGATGGGAAGGGAAGGGAGGGGAAGGGAAGGAACAGGGGAAGGGAAGGGAAGGG".to_string(),
            "AAGGGAAGGGAAGGGAAGGGAAGGGCAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGG".to_string(),
            "AAGGGAAGGGAAGGGATGGGAAGGGCAGGGAAGGGAAGGGAAGGAAAAGGGAAGGGAAGGGAAGGG".to_string(),
            "AAGGGAAGAATGGGAAGGGAAGGGCAGGGAAGGGAAGGGGAAGGGAAGGGAAGGGAAGGGAAGGG".to_string(),
            "AAGGGAAGGGAAGGGATGGGAAGGGCAGGGAAGGGAAGGGAAGGAAAAGGGAAGGGAAGGGAAGGG".to_string(),
            "AAGGGAAGAATGGGAAGGGAAGGGCAGGGAAGGGAAGGGGAAGGGAAGGGAAGGGAAGGGAAGGG".to_string(),
            "AAGGGGGGAAGGGACGGGAAGGGCAGGGAAGGGAAGGGAAGGAAAAAGGGAAGGGAAGGGAAGGG".to_string(),
        ];
        // shuffle the insertions in a random order
        use rand::seq::SliceRandom;
        let mut rng = rand::rng();
        insertions.shuffle(&mut rng);
        let splitseqs = split(
            &insertions,
            &crate::repeats::RepeatInterval {
                chrom: "chr7".to_string(),
                start: 154654404,
                end: 154654432,
                created: None,
            },
            false,
            0.1,
            3,
        );
        let mut hap1 = splitseqs.hap1;
        let mut hap2 = splitseqs.hap2.unwrap();
        // assert that either hap1 or hap2 is equal to the expected_haplotype
        // but first sort the haplotypes, since the order is not guaranteed
        hap1.sort();
        hap2.sort();
        expected_haplotype.sort();
        assert!(hap1 == expected_haplotype || hap2 == expected_haplotype);
    }
}
