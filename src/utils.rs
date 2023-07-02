use flate2::read;
use rust_htslib::bam;
use rust_htslib::bam::record::Aux;
use rust_htslib::faidx;
use std::ffi::OsStr;
use std::fs::File;
use std::io::Read;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Read normal or compressed files seamlessly
/// Uses the presence of a `.gz` extension to decide
pub fn reader(filename: &str) -> Box<dyn BufRead> {
    let path = Path::new(filename);
    let file = match File::open(path) {
        Err(why) => panic!("couldn't open {}: {}", path.display(), why),
        Ok(file) => file,
    };

    if path.extension() == Some(OsStr::new("gz")) {
        Box::new(BufReader::with_capacity(
            128 * 1024,
            read::GzDecoder::new(file),
        ))
    } else {
        Box::new(BufReader::with_capacity(128 * 1024, file))
    }
}

/// parse a region string
pub fn process_region(reg: String) -> Result<(String, u32, u32), Box<dyn std::error::Error>> {
    let chrom = reg.split(':').collect::<Vec<&str>>()[0];
    let interval = reg.split(':').collect::<Vec<&str>>()[1];
    let start: u32 = interval.split('-').collect::<Vec<&str>>()[0]
        .parse()
        .unwrap();
    let end: u32 = interval.split('-').collect::<Vec<&str>>()[1]
        .parse()
        .unwrap();
    assert!(
        end - start > 0,
        r#"Invalid region: begin has to be smaller than end."#
    );
    Ok((chrom.to_string(), start, end))
}

pub fn get_phase(record: &bam::Record) -> u8 {
    match record.aux(b"HP") {
        Ok(value) => {
            if let Aux::U8(v) = value {
                v
            } else {
                panic!("Unexpected type of Aux {value:?}")
            }
        }
        Err(_e) => 0,
    }
}

pub fn write_vcf_header(fasta: &str, bam: &str, sample: Option<String>) {
    println!(r#"##fileformat=VCFv4.2"#);
    // get absolute path to fasta file
    let path = std::fs::canonicalize(fasta)
        .unwrap_or_else(|err| panic!("Failed getting absolute path to fasta: {err}"));
    println!(
        r#"##reference={}"#,
        path.to_str().expect("Failed converting path to string")
    );
    // get the version of this crate
    let version = env!("CARGO_PKG_VERSION");
    println!(r#"##source=STRdust v{}"#, version);
    // call faidx to make sure the fasta index exists, we'll need this anyway when genotyping
    let _ =
        faidx::Reader::from_path(fasta).unwrap_or_else(|err| panic!("Failed opening fasta: {err}"));

    let mut fai_file = std::fs::File::open(format!("{fasta}.fai")).expect("Can't open file");
    // parse the fasta index file
    let mut buf = String::new();
    fai_file
        .read_to_string(&mut buf)
        .expect("Can't read fai file");
    for contig in buf.lines() {
        let mut contig = contig.split_whitespace();
        let name = contig.next().unwrap();
        let length = contig.next().unwrap().parse::<usize>().unwrap();
        println!(r#"##contig=<ID={},length={}>"#, name, length);
    }
    println!(
        r#"##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the repeat interval">"#
    );
    println!(
        r#"##INFO=<ID=RL,Number=2,Type=Integer,Description="Repeat length of the two alleles">"#
    );
    println!(
        r#"##INFO=<ID=SUPP,Number=2,Type=Integer,Description="Number of reads supporting the two alleles">"#
    );
    println!(
        r#"##INFO=<ID=STDEV,Number=2,Type=Integer,Description="Standard deviation of the repeat length">"#
    );
    println!(
        r#"##INFO=<ID=SEQS,Number=1,Type=String,Description="Sequences supporting the two alleles">"#
    );
    println!(r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#);
    println!(r#"##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set identifier">"#);
    let name = match sample {
        Some(name) => name,
        None => {
            // use basename of bam and remove file extension
            let name = std::path::Path::new(&bam)
                .file_stem()
                .unwrap()
                .to_str()
                .unwrap();
            name.to_string()
        }
    };
    println!("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{name}",);
}

#[cfg(test)]
mod tests {
    use rust_htslib::bam::Read;

    use super::*;

    #[test]
    fn test_get_phase() {
        let mut bam =
            bam::Reader::from_path("test_data/small-test-phased.bam").expect("Failed opening bam");
        let record = bam
            .records()
            .next()
            .expect("Failed to read first record from bam");
        let phase = get_phase(&record.unwrap());
        assert_eq!(phase, 2);
    }

    #[test]
    fn test_write_vcf_header_from_bam() {
        write_vcf_header(
            "test_data/chr7.fa.gz",
            "test_data/small-test-phased.bam",
            None,
        );
    }

    #[test]
    fn test_write_vcf_header_from_name() {
        write_vcf_header(
            "test_data/chr7.fa.gz",
            "test_data/small-test-phased.bam",
            Some("test_sample".to_string()),
        );
    }
}
