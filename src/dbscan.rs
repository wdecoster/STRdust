//! DBSCAN clustering for the optional `--phasing-strategy dbscan` path.
//!
//! Adapted from the `trout` cohort-outlier tool, which used a custom DBSCAN
//! over a precomputed CSR adjacency for speed (one pairwise distance per pair,
//! no per-range-query allocations). Here it returns full cluster *labels*
//! (not just a noise mask) because phasing needs the read groupings.
//!
//! Distances are squared Euclidean over the feature vectors from
//! [`crate::features`]; `min_samples` is the standard DBSCAN core-point
//! threshold (a point counts itself, so a core point has ≥ `min_samples`
//! members including itself).

/// Cluster label for each point. `None` == noise (DBSCAN outlier).
pub type Label = Option<usize>;

/// Run DBSCAN and return a per-point label vector. Clustered points get a
/// `Some(cluster_id)`; noise points get `None`. Cluster ids are assigned in
/// the order clusters are discovered (0, 1, 2, …).
pub fn cluster(points: &[Vec<f64>], eps: f64, min_samples: usize) -> Vec<Label> {
    let n = points.len();
    if n == 0 {
        return Vec::new();
    }
    let eps_sq = eps * eps;

    // Collect neighbour edges once, then build CSR adjacency from a single
    // allocation pair (avoids growing n separate small Vecs).
    let mut edges: Vec<(u32, u32)> = Vec::new();
    for i in 0..n {
        let pi = points[i].as_slice();
        for (j, pj_vec) in points.iter().enumerate().skip(i + 1) {
            let pj = pj_vec.as_slice();
            // LLVM auto-vectorises this zipped sum on f64 slices.
            let dsq: f64 = pi
                .iter()
                .zip(pj.iter())
                .map(|(a, b)| {
                    let d = a - b;
                    d * d
                })
                .sum();
            if dsq <= eps_sq {
                edges.push((i as u32, j as u32));
            }
        }
    }

    // CSR build: row_starts[i..i+1] indexes col_indices for point i's neighbours.
    let mut row_lens = vec![1usize; n]; // each point includes itself
    for &(i, j) in &edges {
        row_lens[i as usize] += 1;
        row_lens[j as usize] += 1;
    }
    let mut row_starts = vec![0usize; n + 1];
    for i in 0..n {
        row_starts[i + 1] = row_starts[i] + row_lens[i];
    }
    let mut col_indices: Vec<u32> = vec![0u32; row_starts[n]];
    let mut cursor = row_starts.clone();
    for i in 0..n {
        col_indices[cursor[i]] = i as u32; // self-loop
        cursor[i] += 1;
    }
    for &(i, j) in &edges {
        col_indices[cursor[i as usize]] = j;
        cursor[i as usize] += 1;
        col_indices[cursor[j as usize]] = i;
        cursor[j as usize] += 1;
    }
    drop(edges);
    drop(cursor);
    drop(row_lens);

    let neighbours_of = |i: usize| -> &[u32] { &col_indices[row_starts[i]..row_starts[i + 1]] };

    #[derive(Clone, Copy, PartialEq, Eq)]
    enum St {
        Unvisited,
        Noise,
        Clustered(usize),
    }
    let mut status = vec![St::Unvisited; n];
    let mut queue: Vec<usize> = Vec::new();
    let mut next_cluster = 0usize;

    for p in 0..n {
        if status[p] != St::Unvisited {
            continue;
        }
        if neighbours_of(p).len() < min_samples {
            status[p] = St::Noise;
            continue;
        }
        // p is a core point — start a new cluster and expand.
        let cid = next_cluster;
        next_cluster += 1;
        status[p] = St::Clustered(cid);
        queue.clear();
        for &q in neighbours_of(p) {
            let q = q as usize;
            if q != p && matches!(status[q], St::Unvisited | St::Noise) {
                status[q] = St::Clustered(cid);
                queue.push(q);
            }
        }
        let mut idx = 0;
        while idx < queue.len() {
            let q = queue[idx];
            idx += 1;
            // Only core points propagate the cluster.
            if neighbours_of(q).len() >= min_samples {
                for &r in neighbours_of(q) {
                    let r = r as usize;
                    if matches!(status[r], St::Unvisited | St::Noise) {
                        status[r] = St::Clustered(cid);
                        queue.push(r);
                    }
                }
            }
        }
    }

    status
        .into_iter()
        .map(|s| match s {
            St::Clustered(cid) => Some(cid),
            _ => None,
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty() {
        assert!(cluster(&[], 0.5, 3).is_empty());
    }

    #[test]
    fn two_separated_clusters_plus_noise() {
        // two dense 1-D blobs near 0.0 and 10.0, one far-away noise point
        let mut pts: Vec<Vec<f64>> = Vec::new();
        for i in 0..6 {
            pts.push(vec![0.0 + i as f64 * 0.01]);
        }
        for i in 0..6 {
            pts.push(vec![10.0 + i as f64 * 0.01]);
        }
        pts.push(vec![100.0]); // noise
        let labels = cluster(&pts, 0.5, 3);
        // first 6 share a label, next 6 share another, last is noise
        let c0 = labels[0];
        let c1 = labels[6];
        assert!(c0.is_some() && c1.is_some() && c0 != c1);
        assert!(labels[..6].iter().all(|l| *l == c0));
        assert!(labels[6..12].iter().all(|l| *l == c1));
        assert_eq!(labels[12], None);
    }

    #[test]
    fn all_noise_when_sparse() {
        let pts: Vec<Vec<f64>> = (0..10).map(|i| vec![(i * 100) as f64]).collect();
        let labels = cluster(&pts, 1.0, 3);
        assert!(labels.iter().all(|l| l.is_none()));
    }
}
