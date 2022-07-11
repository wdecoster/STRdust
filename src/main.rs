use rust_htslib::{bam, bam::Read};




fn main() {
    println!("Hello, world!");
}

fn parse_insertions_from_bam(bampath) {
    let bam = bam::Reader::from_path(bampath).unwrap();
    for r in bam.records() {
        let record = r.unwrap();
        if record.is_reverse() {
            out.write(&record).unwrap();
        }
    }
}