use std::fmt;

use bio::alignment::pairwise::Scoring;
use bio::alignment::poa::*;
use log::debug;
use rand::seq::SliceRandom;

#[derive(Clone)]
pub struct Consensus {
    pub seq: Option<String>,
    pub support: usize,
    pub std_dev: usize,
    pub score: i32,
}

impl Default for Consensus {
    fn default() -> Consensus {
        Consensus {
            seq: None,
            support: 0,
            std_dev: 0,
            score: -1,
        }
    }
}

impl Consensus {
    pub fn format_lengths(&self, start: u32, end: u32) -> (String, String, String, String, String) {
        match &self.seq {
            Some(seq) => (
                // length of the consensus sequence minus the length of the repeat sequence
                (seq.len() as i32 - ((end - start) as i32)).to_string(),
                seq.clone(),
                self.support.to_string(),
                self.std_dev.to_string(),
                self.score.to_string(),
            ),
            None => (
                ".".to_string(),
                ".".to_string(),
                self.support.to_string(),
                ".".to_string(),
                ".".to_string(),
            ),
        }
    }
}
//              Some(seq) => {
//                     // length of the consensus sequence minus the length of the repeat sequence
//             (
//                 consensus.seq.len() as i32 - ((end - start) as i32)).to_string(),
//             consensus.seq.clone(),
//             consensus.support.to_string(),
//             consensus.std_dev.to_string()
//         )
//              };
//             None => (
//                 ".".to_string(),
//                 ".".to_string(),
//                 self.support.to_string(),
//                 ".".to_string(),
//             ),
//     }
// }

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

pub fn consensus(seqs: &[String], support: usize) -> Consensus {
    if seqs.is_empty() {
        return Consensus {
            seq: None,
            support: 0,
            std_dev: 0,
            score: -1,
        };
    }
    let num_reads_ = seqs.len();
    let (seqs, std_dev) = remove_outliers(seqs);
    let num_reads = seqs.len();
    debug!(
        "Kept {}/{} reads after dropping outliers",
        num_reads, num_reads_
    );
    if num_reads < support {
        Consensus {
            seq: None,
            support: num_reads,
            std_dev,
            score: -1,
        }
    } else {
        // if there are more than 20 reads, downsample to 20 before taking the consensus
        // for performance and memory reasons
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
        let score = aligner.global_banded(&consensus, 20).alignment().score;

        Consensus {
            seq: Some(std::str::from_utf8(&consensus).unwrap().to_string()),
            support: num_reads,
            std_dev,
            score,
        }
    }
}

fn remove_outliers(seqs: &[String]) -> (Vec<&String>, usize) {
    // remove sequences that are shorter or longer than two standard deviations from the mean
    // except if the stdev is small
    let lengths = seqs.iter().map(|x| x.len()).collect::<Vec<usize>>();
    debug!("lengths: {:?}", lengths);

    let mean = lengths.iter().sum::<usize>() / lengths.len();
    let variance = lengths
        .iter()
        .map(|x| (*x as isize - mean as isize).pow(2) as usize)
        .sum::<usize>()
        / lengths.len();
    // note that the casting to usize will floor the std_dev to the integer below, rather than properly rounding it to the nearest integer
    // so this keeps the std_dev smaller than it really is, but not by a lot
    // the places where this matter are probably negligible
    let std_dev = (variance as f64).sqrt() as usize;
    debug!("mean: {}, std_dev: {}", mean, std_dev);
    if std_dev < 5 {
        debug!("std_dev < 5, not removing any outliers");
        (seqs.iter().collect::<Vec<&String>>(), std_dev)
    } else {
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
        let cons = consensus(&seqs, 10);
        println!("Consensus: {}", cons.seq.unwrap());
        println!("Num reads: {}", cons.support);
        println!("Std dev: {}", cons.std_dev);
        println!("Consensus score: {}", cons.score);
    }
}
