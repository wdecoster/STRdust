# STRdust

Tandem repeat genotyper for long reads, returning the repeat length and sequence in VCF format.

## Usage

```text
    STRdust [OPTIONS] <FASTA> <BAM>

ARGS:
    <FASTA>    reference genome used for alignment
    <BAM>      bam/cram file to call STRs in (local path or URL)

SPECIFY ONE OF:
    -r, --region <REGION>              region string to genotype expansion in
    -R, --region-file <REGION_FILE>    Bed file with region(s) to genotype expansion(s) in
        --pathogenic                   Genotype the pathogenic STRs from STRchive

OPTIONS:
    -m, --minlen <MINLEN>              minimal length of insertion/deletion operation [default: 5]
    -s, --support <SUPPORT>            minimal number of supporting reads per haplotype [default: 3]
    -t, --threads <THREADS>            Number of parallel threads to use [default: 1]
        --sample <SAMPLE>              Sample name to use in VCF header, if not provided, the bam
                                       file name is used
        --somatic                      Print information on somatic variability
        --unphased                     Reads are not phased, will use hierarchical clustering to
                                       phase expansions
        --find-outliers                Identify poorly supported outlier expansions (only with
                                       --unphased)
        --haploid <HAPLOID>            comma-separated list of haploid (sex) chromosomes
    -h, --help                         Print help information
    -V, --version                      Print version information
```

## Installation

Preferably, for most users, download a ready-to-use binary for your system to add directory on your $PATH from the [releases](https://github.com/wdecoster/STRdust/releases).  
You may have to change the file permissions to execute it with `chmod +x STRdust`

## CITATION

If you use this tool, please consider citing our publication.
