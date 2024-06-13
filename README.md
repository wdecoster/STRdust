# STRdust

Tandem repeat genotyper for long reads, returning the repeat length and sequence in VCF format.

## Usage

### Installation

Preferably, for most users, download a ready-to-use binary for your system to add directory on your $PATH from the [releases](https://github.com/wdecoster/STRdust/releases).  
You may have to change the file permissions to execute it with `chmod +x STRdust`

Alternatively, you can install the tool using cargo:

```bash
git clone https://github.com/wdecoster/STRdust.git
cd STRdust
cargo build --release
```

### Quick start examples

```bash
STRdust -r chr7:154654404-154654432 reference.fa sample.cram > sample.vcf
STRdust --pathogenic reference.fa sample.cram | bgzip > sample.vcf.gz
STRdust -R targets.bed --haploid chrX,chrY reference.fa male_sample.cram | bgzip > repeats.vcf.gz
```

The 'test_data' directory contains a small example dataset to test the tool:

```bash
STRdust -r chr7:154654404-154654432 test_data/chr7.fa.gz test_data/small-test-phased.bam > small-test-phased.vcf
```

### All arguments

```text
    STRdust [OPTIONS] <FASTA> <BAM>

ARGS:
    <FASTA>    reference genome used for alignment
    <BAM>      bam/cram file to call STRs in (local path or URL)

SPECIFY ONE OF:
    -r, --region <REGION>              region string to genotype expansion in (format: chr:start-end)
    -R, --region-file <REGION_FILE>    Bed file with region(s) to genotype expansion(s) in
        --pathogenic                   Genotype the pathogenic STRs from STRchive

OPTIONS:
    -m, --minlen <MINLEN>              minimal length of insertion/deletion operation [default: 5]
    -s, --support <SUPPORT>            minimal number of supporting reads per haplotype [default: 3]
    -t, --threads <THREADS>            Number of parallel threads to use [default: 1]
        --sample <SAMPLE>              Sample name to use in VCF header, if not provided, the bam file name is used
        --somatic                      Print information on somatic variability
        --unphased                     Reads are not phased, will use hierarchical clustering to phase expansions
        --find-outliers                Identify poorly supported outlier expansions (only with --unphased)
        --haploid <HAPLOID>            comma-separated list of haploid (sex) chromosomes
    -h, --help                         Print help information
    -V, --version                      Print version information
```

## Output format

STRdust produces a VCF file per sample. The consensus sequence is in the ALT field, with sequences from each read in the SEQS INFO field (when running with --somatic).
The FRB FORMAT field is the total repeat length, of the two alleles, in nucleotides. The RB field is the difference between the indidiual allele lengths and the reference length.
The SC FORMAT field is a measure of accuracy of the consensus sequence compared to the overlap graph from the individual reads, which could be influenced by the presence of sequencing errors or somatic variation.

Example output:

```text
(header cropped for brevity)
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the repeat interval">
##INFO=<ID=STDEV,Number=2,Type=Integer,Description="Standard deviation of the repeat length">
##INFO=<ID=SEQS,Number=1,Type=String,Description="Sequences supporting the two alleles">
##INFO=<ID=OUTLIERS,Number=1,Type=String,Description="Outlier sequences much longer than the alleles">
##INFO=<ID=CLUSTERFAILURE,Number=0,Type=Flag,Description="If unphased input failed to cluster in two haplotype">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=RB,Number=2,Type=Integer,Description="Repeat length of the two alleles in bases relative to reference">
##FORMAT=<ID=FRB,Number=2,Type=Integer,Description="Full repeat length of the two alleles in bases">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set identifier">
##FORMAT=<ID=SUP,Number=2,Type=Integer,Description="Read support per allele">
##FORMAT=<ID=SC,Number=2,Type=Integer,Description="Consensus score per allele">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  HG00271.hg38
chr1    1435798 .       TGGCGCGGAGCGGCGCGGAGCG  GCTGGCGCGGAGCGGCGCGGA,GCGGGCGCGCGCAGGA  .       .       END=1435818;STDEV=1,2     GT:RB:FRB:SUP:SC        1|2:1,-4:21,16:18,6:63,41
chr1    57367044        .       AAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAATAAAT     AAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAATAAA,AAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAATAAA        .       .       END=57367125;STDEV=3,0 GT:RB:FRB:SUP:SC 1|2:-9,-34:72,47:17,12:216,141
```

## CITATION

If you use this tool, please consider citing our [publication](https://www.medrxiv.org/content/10.1101/2024.03.06.24303700v1).
