// write a function that takes a long string like CAGCAGCAGCAGCGGCGGCGGCAGCAGCAG and converts it to a condensed representation like (CAG)4(CGG)3(CAG)3
// the function should identify repeated parts and replace them with (repeat)number

#[allow(dead_code)]
fn create_motif(_seq: &str) -> String {
    unimplemented!()
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
}
