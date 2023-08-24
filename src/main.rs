#![allow(non_snake_case)]
use clap::AppSettings::DeriveDisplayOrder;
use clap::Parser;
use log::info;
use std::path::PathBuf;

pub mod call;
pub mod consensus;
pub mod motif;
pub mod phase_insertions;
pub mod repeats;
pub mod utils;
pub mod vcf;

// The arguments end up in the Cli struct
#[derive(Parser, Debug)]
#[structopt(global_settings=&[DeriveDisplayOrder])]
#[clap(author, version, about="Tool to genotype STRs from long reads", long_about = None)]
pub struct Cli {
    /// reference genome
    #[clap(parse(from_os_str), validator=is_file)]
    fasta: PathBuf,

    /// bam file to call STRs in
    #[clap(parse(from_os_str), validator=is_file)]
    bam: PathBuf,

    /// region string to genotype expansion in
    #[clap(short, long, value_parser)]
    region: Option<String>,

    /// Bed file with region(s) to genotype expansion(s) in
    #[clap(short = 'R', long, value_parser, validator=is_file)]
    region_file: Option<PathBuf>,

    /// minimal length of insertion/deletion operation
    #[clap(short, long, value_parser, default_value_t = 5)]
    minlen: usize,

    /// minimal number of supporting reads per haplotype
    #[clap(short, long, value_parser, default_value_t = 3)]
    support: usize,

    /// Number of parallel threads to use
    #[clap(short, long, value_parser, default_value_t = 8)]
    threads: usize,

    /// Sample name to use in VCF header, if not provided, the bam file name is used
    #[clap(long, value_parser)]
    sample: Option<String>,

    /// Print information on somatic variability
    #[clap(long, value_parser, default_value_t = false)]
    somatic: bool,

    /// Reads are not phased
    #[clap(long, value_parser, default_value_t = false)]
    unphased: bool,
}

fn is_file(pathname: &str) -> Result<(), String> {
    let path = PathBuf::from(pathname);
    if path.is_file() {
        Ok(())
    } else {
        Err(format!("Input file {} is invalid", path.display()))
    }
}

fn main() {
    env_logger::init();
    let args = Cli::parse();
    info!("Collected arguments");
    call::genotype_repeats(args);
}

#[cfg(test)]
#[ctor::ctor]
fn init() {
    env_logger::init();
}

#[test]
fn verify_app() {
    use clap::CommandFactory;
    Cli::command().debug_assert()
}
