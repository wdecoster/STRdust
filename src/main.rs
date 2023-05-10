#![allow(non_snake_case)]
use clap::AppSettings::DeriveDisplayOrder;
use clap::Parser;
use log::info;
use std::path::PathBuf;

pub mod call;
pub mod consensus;
pub mod utils;

// The arguments end up in the Cli struct
#[derive(Parser, Debug)]
#[structopt(global_settings=&[DeriveDisplayOrder])]
#[clap(author, version, about="Tool to genotype STRs from long reads", long_about = None)]
struct Cli {
    /// bam file to call STRs in
    #[clap(parse(from_os_str), validator=is_file)]
    bam: PathBuf,

    /// reference genome
    #[clap(parse(from_os_str), validator=is_file)]
    fasta: PathBuf,

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

    /// Print information on somatic variability
    #[clap(long, value_parser, default_value_t = false)]
    somatic: bool,
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
    call::genotype_repeats(
        args.bam,
        args.fasta,
        args.region,
        args.region_file,
        args.minlen,
        args.support,
        args.threads,
        args.somatic,
    );
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
