// write a function that takes a long string like CAGCAGCAGCAGCGGCGGCGGCAGCAGCAG and converts it to a condensed representation like (CAG)4(CGG)3(CAG)3
// the function should identify repeated parts and replace them with (repeat)number

fn create_motif(seq: &str) -> String {
    unimplemented!()
}

pub fn create_mc(motifs_string: &str, repeats: &str) -> String {
     let motifs: Vec<String> = motifs_string.split(',').map(|s| s.to_string()).collect();

     let mut motif_counts = vec![0; motifs.len()];
     for motif in motifs.iter() {
         let mut start = 0;
         while let Some(index) = repeats[start..].find(motif) {
             motif_counts[motifs.iter().position(|x| x == motif).unwrap()] += 1;
             start += index + motif.len();
         }
     }
     motif_counts.iter().map(|count| count.to_string()).collect::<Vec<String>>().join("_")
 }

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[ignore]
    fn test_create_motif() {
        assert_eq!(
            create_motif("CAGCAGCAGCAGCGGCGGCGGCAGCAGCAG"),
            "(CAG)4(CGG)3(CAG)3"
        );
    }
    #[test]
    fn test_create_mc() {
        let motifs = "ATG,CGT,GCA";
        let repeat = "ATGCGTATGCGTAGCGT";
        println!("{:?}", create_mc(&motifs, &repeat));
        assert_eq!(
            create_mc(&motifs, repeat), "2_3_0"
        );
    }
    #[test]
    fn test_create_mc_with_n() {
        let motifs = "NCG";
        let repeat = "ACGACGACGACG";

        assert_eq!(
            create_mc(&motifs, repeat), "4"
        );
    }
}
