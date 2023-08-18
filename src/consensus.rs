use bio::alignment::pairwise::Scoring;
use bio::alignment::poa::*;
use log::{debug, log_enabled, Level};
use rand::seq::SliceRandom;
use std::error::Error;
use std::fmt;

#[derive(Debug, Clone)]
pub struct ConsensusError {
    pub support: usize,
}

impl ConsensusError {
    fn new(support: usize) -> Self {
        ConsensusError { support }
    }
}

impl fmt::Display for ConsensusError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "Not enough reads to generate consensus: {}",
            self.support
        )
    }
}

impl Error for ConsensusError {
    fn description(&self) -> &str {
        "Not enough reads to generate consensus"
    }
}

pub fn consensus(
    seqs: &[String],
    support: usize,
) -> Result<(String, usize, usize), ConsensusError> {
    if seqs.is_empty() {
        return Err(ConsensusError::new(0));
    }
    let num_reads_ = seqs.len();
    let (seqs, std_dev) = remove_outliers(seqs);
    let num_reads = seqs.len();
    debug!(
        "Kept {}/{} reads after dropping outliers",
        num_reads, num_reads_
    );
    if num_reads < support {
        Err(ConsensusError::new(num_reads))
    } else {
        // if there are more than 20 reads, downsample to 20 before taking the consensus
        let seqs = if num_reads > 20 {
            debug!("Too many reads, downsampling to 20");
            seqs.choose_multiple(&mut rand::thread_rng(), 20)
                .cloned()
                .collect::<Vec<&String>>()
        } else {
            seqs
        };
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
            aligner.global_banded(seq, 20).add_to_graph();
        }
        let consensus = aligner.consensus();

        Ok((
            std::str::from_utf8(&consensus).unwrap().to_string(),
            num_reads,
            std_dev,
        ))
    }
}

fn remove_outliers(seqs: &[String]) -> (Vec<&String>, usize) {
    // remove sequences that are shorter or longer than two standard deviations from the mean
    // except if the stdev is small
    let mut lengths = seqs.iter().map(|x| x.len()).collect::<Vec<usize>>();
    if log_enabled!(Level::Debug) {
        lengths.sort();
        debug!("lengths: {:?}", lengths);
    }
    let mean = lengths.iter().sum::<usize>() / lengths.len();
    let variance = lengths
        .iter()
        .map(|x| (*x as isize - mean as isize).pow(2) as usize)
        .sum::<usize>()
        / lengths.len();
    let std_dev = (variance as f64).sqrt() as usize;
    debug!("mean: {}, std_dev: {}", mean, std_dev);
    if std_dev < 5 {
        debug!("std_dev < 5, not removing any outliers");
        return (seqs.iter().collect::<Vec<&String>>(), std_dev);
    }
    // avoid underflowing usize
    let min_val = mean.saturating_sub(2 * std_dev);
    let max_val = mean + 2 * std_dev;
    debug!("Removing outliers outside of [{},{}]", min_val, max_val);
    let filtered_seqs = seqs
        .iter()
        .zip(lengths.iter())
        .filter(|(_, &len)| len > min_val && len < max_val)
        .map(|(seq, _)| seq)
        .collect::<Vec<&String>>();
    (filtered_seqs, std_dev)
}

#[cfg(test)]
mod tests {
    use super::*;
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
        let (consensus_seq, num_reads, std_dev) = consensus(&seqs, 10).unwrap();
        println!("Consensus: {}", consensus_seq);
        println!("Num reads: {}", num_reads);
        println!("Std dev: {}", std_dev);
    }
}
