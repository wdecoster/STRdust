use log::debug;
use rust_spoa::poa_consensus;

pub fn consensus(seqs: &[String], support: usize) -> Option<(String, usize, usize)> {
    if seqs.is_empty() {
        return None;
    }
    let (seqs, std_dev) = remove_outliers(seqs);
    let num_reads = seqs.len();
    debug!("{} reads after removing outliers", num_reads);
    if num_reads < support {
        None
    } else {
        let consensus_max_length = seqs.iter().map(|x| x.len()).max().unwrap_or(0);
        debug!("{} is the longest read", consensus_max_length);
        let mut seqs_bytes = vec![];
        for seq in seqs.iter() {
            seqs_bytes.push(format!("{seq}\0").bytes().collect::<Vec<u8>>());
        }
        // I empirically determined the following parameters to be suitable,
        // but further testing on other repeats would be good
        // mainly have to make sure the consensus does not get longer than the individual insertions
        let consensus = poa_consensus(
            &seqs_bytes,
            consensus_max_length,
            1,   // 0 = local, 1 = global, 2 = gapped
            2,   // match_score,
            -4,  // mismatch_score,
            -12, // gap_open,
            -6,  // gap_extend,
        );
        Some((
            std::str::from_utf8(&consensus).unwrap().to_string(),
            num_reads,
            std_dev,
        ))
    }
}

fn remove_outliers(seqs: &[String]) -> (Vec<&String>, usize) {
    // remove sequences that are shorter or longer than two standard deviations from the mean
    let lengths = seqs.iter().map(|x| x.len()).collect::<Vec<usize>>();
    debug!("lengths: {:?}", lengths);
    let mean = lengths.iter().sum::<usize>() / lengths.len();
    debug!("mean: {}", mean);
    let variance = lengths
        .iter()
        .map(|x| (*x as isize - mean as isize).pow(2) as usize)
        .sum::<usize>()
        / lengths.len();
    let std_dev = (variance as f64).sqrt() as usize;
    debug!("std_dev: {}", std_dev);
    // avoid underflowing usize
    let min_val = mean.saturating_sub(2 * std_dev);
    debug!("min_val: {}", min_val);
    (
        seqs.iter()
            .zip(lengths.iter())
            .filter(|(_, &len)| len > min_val && len < mean + 2 * std_dev)
            .map(|(seq, _)| seq)
            .collect::<Vec<&String>>(),
        std_dev,
    )
}
