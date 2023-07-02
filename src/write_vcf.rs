use distance::levenshtein;
use log::debug;

pub fn write_vcf(
    mut consenses: Vec<Option<(String, usize, usize)>>,
    repeat_ref_sequence: String,
    somatic: bool,
    all_insertions: Vec<String>,
    chrom: String,
    start: u32,
    end: u32,
) -> String {
    // since I use .pop() to format the two consensus sequences, the order is reversed
    let (length2, alt2, support2, std_dev2) = format_lengths(consenses.pop().unwrap(), start, end);
    let (length1, alt1, support1, std_dev1) = format_lengths(consenses.pop().unwrap(), start, end);
    debug!("Genotyping {chrom}:{start}-{end}:{repeat_ref_sequence} with {alt1} and {alt2}");
    let allele1 = if alt1 == "." {
        "."
    } else if levenshtein(&alt1, &repeat_ref_sequence) < 5 {
        "0"
    } else {
        "1"
    };

    let allele2 = if alt2 == "." {
        "."
    } else if levenshtein(&alt2, &repeat_ref_sequence) < 5 {
        "0"
    } else if alt2 == alt1 {
        "1"
    } else {
        "2"
    };

    let alts = match (allele1, allele2) {
        ("1", "0") | ("1", ".") | ("1", "1") => alt1, // if both alleles are the same, only report one
        ("0", "1") | (".", "1") => alt2,
        ("1", "2") => alt1 + "," + &alt2,
        _ => ".".to_string(), // includes ./. and 0/0
    };

    let somatic_info_field = if somatic {
        format!(";SEQS={}", all_insertions.join("|"))
    } else {
        "".to_string()
    };

    format!(
            "{chrom}\t{start}\t.\t{repeat_ref_sequence}\t{alts}\t.\t.\t\
            END={end};RL={length1}|{length2};SUPP={support1}|{support2};STDEV={std_dev1}|{std_dev2}{somatic_info_field}\
            \tGT\t{allele1}|{allele2}",
        )
}

fn format_lengths(
    consensus: Option<(String, usize, usize)>,
    start: u32,
    end: u32,
) -> (String, String, String, String) {
    match consensus {
        Some((ref s, sup, std_dev)) => (
            // length of the consensus sequence minus the length of the repeat sequence
            (s.len() as i32 - ((end - start) as i32)).to_string(),
            s.clone(),
            sup.to_string(),
            std_dev.to_string(),
        ),
        None => (
            ".".to_string(),
            ".".to_string(),
            ".".to_string(),
            ".".to_string(),
        ),
    }
}

// pub fn empty_record(
//     reason: &str,
//     chrom: String,
//     start: u32,
//     end: u32,
//     repeat_ref_sequence: String,
// ) -> String {
//     eprintln!("Cannot genotype repeat at {chrom}:{start}-{end} because {reason}");
//     format!("{chrom}\t{start}\t.\t{ref}\t.,.\t.\t.\tEND={end};RL=.|.;SUPP=.|.;STDEV=.|.;\tGT\t.|.", ref = repeat_ref_sequence,)
// }
