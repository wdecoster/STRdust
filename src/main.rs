use clap::AppSettings::DeriveDisplayOrder;
use clap::{Parser, Subcommand};
use log::info;
use std::path::PathBuf;

pub mod call;
pub mod utils;

// The arguments end up in the Cli struct
#[derive(Parser, Debug)]
#[structopt(global_settings=&[DeriveDisplayOrder])]
#[clap(author, version, about="Tool to genotype STRs from long reads", long_about = None)]
struct Cli {
    #[clap(subcommand)]
    command: Commands,
}
// Every subcommand is a variation of the Commands Enum, and has its arguments defined below
#[derive(Debug, Subcommand)]
enum Commands {
    /// Call lengths
    #[clap(arg_required_else_help = true)]
    Call {
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
        minlen: u32,

        /// minimal number of supporting reads
        #[clap(short, long, value_parser, default_value_t = 3)]
        support: usize,

        /// Number of parallel threads to use
        #[clap(short, long, value_parser, default_value_t = 8)]
        threads: usize,
    },
    /// Combine lengths from multiple bams to a TSV
    Combine {},
    /// Search for regions potentially containing a polymorphic repeat
    Scan {},
    /// Find outliers from TSV
    Outlier {},
    /// Lookup genotypes and display
    Query {},
    /// Test for association of repeat length by comparing two cohorts
    Association {},
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
    match args.command {
        Commands::Call {
            bam,
            fasta,
            region,
            region_file,
            minlen,
            support,
            threads,
        } => call::genotype_repeats(bam, fasta, region, region_file, minlen, support, threads),
        Commands::Combine {} => {
            unimplemented!();
        }
        Commands::Scan {} => {
            unimplemented!();
        }
        Commands::Outlier {} => {
            unimplemented!();
        }
        Commands::Query {} => {
            unimplemented!();
        }
        Commands::Association {} => {
            unimplemented!();
        }
    }
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
