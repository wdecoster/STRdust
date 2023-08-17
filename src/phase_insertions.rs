use kodama::{linkage, Method};
use levenshtein::levenshtein;
use log::{debug, error, log_enabled, Level};
use std::collections::HashMap;

pub fn split(insertions: &Vec<String>) -> (Option<Vec<String>>, Option<Vec<String>>) {
    // the insertions are from an unphased experiment
    // and should be split in one (if homozygous) or two haplotypes
    // this is based on the length of the insertion
    // as well as the sequence composition
    // the sequence composition has to be agnostic to the motif length
    // and will use the levenshtein distance
    // this is inspired by the TRGT paper

    // Create a condensed (upper triangle) distance matrix
    let mut condensed = vec![];
    for row in 0..insertions.len() - 1 {
        for col in row + 1..insertions.len() {
            condensed.push(levenshtein(&insertions[row], &insertions[col]) as f32);
        }
    }
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
    // clusters have to represent at least 20% of the reads
    let min_cluster_size = (insertions.len() as f32 / 5.0) as usize;
    debug!("Minimum cluster size: {}", min_cluster_size);

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
                "Node {cluster} with dissimilarity {} and children {} [{seq1}] and {} [{seq2}]",
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
    debug!("Roots for this tree: {:?}", roots);

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
                debug!("Adding cluster {} to candidate haplotype clusters", cluster);
            }
            large_cluster_seen.push(*cluster);
        }
    }
    match haplotype_clusters.len() {
        0 => {
            error!("No haplotype clusters found");
            panic!();
        }
        1 => {
            debug!("Only one haplotype cluster found");
            (
                Some(find_cluster_members(
                    &haplotype_clusters[0],
                    &cluster_to_subclusters,
                    insertions,
                )),
                None,
            )
        }
        2 => {
            debug!("Found two haplotype clusters");
            (
                Some(find_cluster_members(
                    &haplotype_clusters[0],
                    &cluster_to_subclusters,
                    insertions,
                )),
                Some(find_cluster_members(
                    &haplotype_clusters[1],
                    &cluster_to_subclusters,
                    insertions,
                )),
            )
        }
        _ => {
            error!("Found more than two haplotype clusters");
            panic!();
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
            roots.push(*child1);
            roots.extend(find_roots(
                *child1,
                cluster_to_subclusters,
                clusters_to_size,
                clusters_to_dissimilarity,
                min_cluster_size,
            ));
        } else {
            roots.push(*child2);
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

fn find_cluster_members(
    cluster: &usize,
    cluster_to_subclusters: &HashMap<usize, (usize, usize)>,
    insertions: &Vec<String>,
) -> Vec<String> {
    let mut search_cluster = vec![*cluster];
    let mut insertions_of_this_cluster = vec![];
    loop {
        let (cluster1, cluster2) = cluster_to_subclusters
            .get(&search_cluster.pop().unwrap())
            .unwrap_or_else(|| panic!("Cluster not in hashmap"));
        for cl in &[cluster1, cluster2] {
            if cl < &&insertions.len() {
                insertions_of_this_cluster.push(insertions[**cl].clone());
            } else {
                search_cluster.push(**cl);
            }
        }
        if search_cluster.is_empty() {
            break;
        }
    }

    insertions_of_this_cluster
}

mod tests {
    #[allow(unused_imports)]
    use super::*;

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
        let mut rng = rand::thread_rng();
        insertions.shuffle(&mut rng);
        let (hap1, hap2) = split(&insertions);
        let hap1 = hap1.unwrap();
        let hap2 = hap2.unwrap();
        assert!(hap1.len() == hap2.len());
        // check that all sequences in hap1 are the same length
        assert!(hap1.iter().all(|x| x.len() == hap1[0].len()));
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
        let mut rng = rand::thread_rng();
        insertions.shuffle(&mut rng);
        let (hap1, hap2) = split(&insertions);
        let mut hap1 = hap1.unwrap();
        let mut hap2 = hap2.unwrap();
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
        let mut rng = rand::thread_rng();
        insertions.shuffle(&mut rng);
        let (hap1, hap2) = split(&insertions);
        let hap1 = hap1.unwrap();
        let hap2 = hap2.unwrap();
        assert!(hap1.len() + hap2.len() == insertions.len());
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
        let mut rng = rand::thread_rng();
        insertions.shuffle(&mut rng);
        let (hap1, hap2) = split(&insertions);
        let hap1 = hap1.unwrap();
        assert!(hap2.is_none());
        println!("hap1: {:?}", hap1);
        assert!(hap1.len() == insertions.len());
    }

    #[test]
    fn test_split_with_outlier() {
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
        let mut rng = rand::thread_rng();
        insertions.shuffle(&mut rng);
        let (hap1, hap2) = split(&insertions);
        let mut hap1 = hap1.unwrap();
        let mut hap2 = hap2.unwrap();
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
    fn test_split_with_two_outliers() {
        // test that the split function identifies two haplotype, based on sequence composition
        // using a AAGGG and AAAAG expansion with some SNV and length noise
        // this is the composition in the RFC1 expansion in CANVAS, with the former being pathogenic
        // additionally, two outlier insertions are added that do not belong to any of the haplotypes
        // these outliers should be discarded
        // additional 'good' reads are added to increase the coverage
        let mut insertions = vec![
            "AAGGGAAGGGAAGGGAAGGGAATGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGG".to_string(),
            "AAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGAAGGGAAGGG".to_string(),
            "AAGGGAAGGGAAGGGAAGGGAAGGGCAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGG".to_string(),
            "AAGGGAAGGGAAGGGATGGGAAGGGCAGGGAAGGGAAGGGAAGGAAAAGGGAAGGGAAGGGAAGGG".to_string(),
            "AAAAGAAAAGAAAAGAAAAGAAAAGAAAGGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAG".to_string(),
            "AAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAGGAAAAGAAAAGAAGAAAAGAAAAGAAAAG".to_string(),
            "AAAAGAAAAGAAAAGAAAAGAAAAAGAAAAGAAAAGAAAAGAAACGAAAAGAAAAGAAAAGAAAAG".to_string(),
            "AAAAGAAAGAAAGAAAAGAAAAAGGGAAAAGAAAAGAAAAGAAACGAAAAGAAAAGAAAAGAAAAAG".to_string(),
            "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT".to_string(),
            "CCCCCCCCCCCCCCCCCCCCCC".to_string(),
        ];
        // the expected_haplotype are the first three insertions
        let mut expected_haplotype = vec![
            "AAGGGAAGGGAAGGGAAGGGAATGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGG".to_string(),
            "AAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGAAGGGAAGGG".to_string(),
            "AAGGGAAGGGAAGGGAAGGGAAGGGCAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGG".to_string(),
            "AAGGGAAGGGAAGGGATGGGAAGGGCAGGGAAGGGAAGGGAAGGAAAAGGGAAGGGAAGGGAAGGG".to_string(),
        ];
        // shuffle the insertions in a random order
        use rand::seq::SliceRandom;
        let mut rng = rand::thread_rng();
        insertions.shuffle(&mut rng);
        let (hap1, hap2) = split(&insertions);
        let mut hap1 = hap1.unwrap();
        let mut hap2 = hap2.unwrap();
        assert!(hap1.len() == hap2.len());
        // assert that either hap1 or hap2 is equal to the expected_haplotype
        // but first sort the haplotypes, since the order is not guaranteed
        hap1.sort();
        hap2.sort();
        expected_haplotype.sort();
        assert!(hap1 == expected_haplotype || hap2 == expected_haplotype);
    }
}
