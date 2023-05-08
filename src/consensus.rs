use rust_spoa::poa_consensus;

pub fn consensus(seqs: &[String], support: usize) -> Option<String> {
    // remove sequences that are shorter or longer than two standard deviations from the mean
    let lengths = seqs.iter().map(|x| x.len()).collect::<Vec<usize>>();
    let mean = lengths.iter().sum::<usize>() / lengths.len();
    let variance = lengths
        .iter()
        .map(|x| (*x as isize - mean as isize).pow(2) as usize)
        .sum::<usize>()
        / lengths.len();
    let std_dev = (variance as f64).sqrt() as usize;
    let seqs = seqs
        .iter()
        .zip(lengths.iter())
        .filter(|(_, &len)| len > mean - 2 * std_dev && len < mean + 2 * std_dev)
        .map(|(seq, _)| seq)
        .collect::<Vec<&String>>();
    if seqs.len() < support {
        None
    } else {
        let consensus_max_length = seqs.iter().map(|x| x.len()).max().unwrap_or(0);
        let mut seqs_bytes = vec![];
        for seq in seqs.iter() {
            seqs_bytes.push(format!("{seq}\0").bytes().collect::<Vec<u8>>());
        }
        let alignment_type = 1; // 0 = local, 1 = global, 2 = gapped
        let match_score = 5;
        let mismatch_score = -4;
        let gap_open = -3;
        let gap_extend = -1;

        let consensus = poa_consensus(
            &seqs_bytes,
            consensus_max_length,
            alignment_type,
            match_score,
            mismatch_score,
            gap_open,
            gap_extend,
        );
        Some(
            std::str::from_utf8(&consensus)
                .unwrap()
                .to_string()
                .to_uppercase(),
        )
    }
}
