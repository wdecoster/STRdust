# STRdust

STRdust is a tandem repeat genotyper for long reads, returning the repeat length and sequence in VCF format. STRdust performs realignment of reads overlapping with a repeat locus to an artificial reference sequence from which the repeat was removed. STRdust was developed while investigating the [pathogenic GOLGA8A repeat expansion](https://www.nature.com/articles/s41588-026-02537-7), and as such was not primarily intended as a general-purpose tandem repeat genotyper. This is reflected by [a benchmark of repeat genotypers](https://www.biorxiv.org/content/10.64898/2026.02.28.708646v1), where STRdust v0.16.0 performs the best of the tested tools on the detection of pathogenic alleles, but less so in the nucleotide-level precision of repeat lengths without expansion.

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
    -r, --region <REGION>              region string to genotype expansion in (format: chr:start-end, 1-based inclusive)
    -R, --region-file <REGION_FILE>    Bed file with region(s) to genotype expansion(s) in (supports .bed and .bed.gz)
        --pathogenic                   Genotype the pathogenic STRs from STRchive

OPTIONS:
    -m, --minlen <MINLEN>              minimal length of insertion/deletion operation [default: 1]
    -s, --support <SUPPORT>            minimal number of supporting reads per haplotype [default: 3]
    -t, --threads <THREADS>            Number of parallel threads to use [default: 1]
        --sample <SAMPLE>              Sample name to use in VCF header, if not provided, the bam file name is used
        --somatic                      Print information on somatic variability
        --unphased                     Reads are not phased, will cluster the reads to phase expansions
        --consensus-reads              Maximum number of reads to use to build the consensus sequence [default: 20]
        --max-number-reads             Max number of reads to extract per locus for genotyping (-1 for all reads) [default: 60]
        --max-locus <MAX_LOCUS>        Maximum locus size to consider; larger intervals are filtered out
        --find-outliers                Identify poorly supported outlier expansions (only with --unphased)
        --min-haplotype-fraction <F>   Minimum fraction of reads for a cluster to be a haplotype (only with --unphased) [default: 0.1]
        --phasing <STRATEGY>           How to split unphased reads into haplotypes: 'ward', 'dbscan' or 'both' (only with --unphased) [default: ward]
        --haploid <HAPLOID>            comma-separated list of haploid (sex) chromosomes
        --alignment-all                Always use full alignment (disable fast reference check via CIGAR)
        --sorted                       Sort output by chrom, start and end
        --debug                        Debug mode
    -h, --help                         Print help information
    -V, --version                      Print version information
```

## Notes

- BED files can be provided in plain text or gzipped format (`.bed` or `.bed.gz`)
- Lowering the number of consensus reads may lead to lesser accurate alternative allele sequences (selecting randomly from the reads), but may greatly improve speed. Note that in the case of somatic length variation, a small number of randomly selected reads may lead to a bias and not be representative of the true repeat length.
- Genotyping known pathogenic repeats with the `--pathogenic` flag will return a VCF with the pathogenic STRs from STRchive, but currently only for the GRCh38 reference.
- For unphased data (`--unphased`), STRdust splits the reads into (at most) two haplotypes before building consensus. Three strategies are available via `--phasing`:
  - `ward` (default): hierarchical (Ward) clustering on a **length-weighted Levenshtein distance**. The distance is dominated by how much two reads differ in length, so reads are grouped primarily by length. Robust for the common case where alleles differ mainly in length, but it tends to fragment a single length-variable expansion across several length bins.
  - `dbscan` (experimental): DBSCAN on **length-invariant k-mer composition feature vectors**, i.e. it groups reads by their *sequence composition* (which motif they are made of) rather than primarily by length. This is the key difference from `ward`: two reads of the same motif cluster together even when their lengths differ a lot, so a length-variable expansion is kept together as one allele. The trade-off is that the reference and expanded alleles must differ in composition for DBSCAN to separate them. The two largest clusters become the haplotypes; remaining clusters and noise reads are reported as `OUTLIERS`, and the total number of clusters is reported in `NCLUSTERS` (so loci with `NCLUSTERS > 2`, i.e. complex/multi-population loci, can be flagged downstream). Its internal parameters (neighbourhood radius and length weight) are hardcoded — see [Tuning and hardcoded parameters](#tuning-and-hardcoded-parameters).
  - `both` (QC mode): report the `ward` call as usual, but additionally run `dbscan` and, when the two disagree by more than 2x on the longer allele, raise a `DISCORDANT_LENGTH` flag and report the DBSCAN allele lengths in `DBSCAN_RB`. This deliberately over-flags (it prioritises sensitivity) and is intended as a triage signal: a worklist of loci worth reviewing where the two orthogonal clustering approaches disagree. The reported genotype is unchanged from `ward`.
- By default, STRdust uses a fast reference check (QUICKREF) to skip full alignment at loci that appear to be homozygous reference. It inspects the CIGAR strings of the first 25 reads spanning a locus, and if at least 5 are found and none show a length difference from the reference of more than 3 bp, the locus is called 0|0 immediately. Loci called this way are marked with a `QUICKREF` flag in the VCF INFO field. This substantially speeds up runs on samples with many reference-like loci. To disable this optimisation and always perform full alignment, use `--alignment-all`.

## Output format

STRdust produces a VCF file per sample. The consensus sequence is in the ALT field, with sequences from each read in the SEQS INFO field (when running with --somatic).
The FRB FORMAT field is the total repeat length, of the two alleles, in nucleotides. The RB field is the difference between the indidiual allele lengths and the reference length.
The MRL FORMAT field is the median read length per allele (relative to the reference). Because it is a median rather than the POA consensus length, it is more robust to a long tail of length-variable reads, and comparing it to RB highlights alleles whose consensus length is being pulled by a skewed read distribution.
The SC FORMAT field is a measure of accuracy of the consensus sequence compared to the overlap graph from the individual reads, which could be influenced by the presence of sequencing errors or somatic variation.

Two locus-level quality flags are raised automatically when running with `--unphased` (no extra options needed):

- `EXPANSION_OUTLIER`: at least 2 reads are more than 2x longer than the longer called allele. This catches a larger expansion that was missed by both alleles because its few reads fell to noise/outliers during clustering, or were dropped as length outliers when building the consensus. The thresholds (2 reads, 2x) are currently fixed; they were chosen to suppress single-read artifacts while still surfacing genuine multi-read expansions.
- `IMPRECISE_LENGTH`: a called allele's reads have a wide length spread (coefficient of variation, std_dev/mean, above 0.2). At such loci the underlying length distribution is continuous or long-tailed, so a single consensus length is not a faithful summary — compare RB and MRL, and treat the call with caution. The 0.2 threshold is currently fixed.

Example output:

```text
(header cropped for brevity)
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the repeat interval">
##INFO=<ID=STDEV,Number=2,Type=Integer,Description="Standard deviation of the repeat length">
##INFO=<ID=SEQS,Number=1,Type=String,Description="Sequences supporting the two alleles">
##INFO=<ID=OUTLIERS,Number=1,Type=String,Description="Outlier sequences much longer than the alleles">
##INFO=<ID=NCLUSTERS,Number=1,Type=Integer,Description="Number of read clusters found by the DBSCAN phasing strategy (>2 indicates a complex multi-population locus)">
##INFO=<ID=CLUSTERFAILURE,Number=0,Type=Flag,Description="If unphased input failed to cluster in two haplotype">
##INFO=<ID=EXPANSION_OUTLIER,Number=0,Type=Flag,Description="At least 2 reads are more than 2x longer than the longer called allele, suggesting a larger expansion missed by both alleles (only with --unphased)">
##INFO=<ID=IMPRECISE_LENGTH,Number=0,Type=Flag,Description="A called allele has a wide read-length spread (coefficient of variation > 0.2); the reported consensus length may not be representative">
##INFO=<ID=DISCORDANT_LENGTH,Number=0,Type=Flag,Description="With --phasing both: the reported Ward call and the DBSCAN call differ by more than 2x on the longer allele; see DBSCAN_RB (QC only, fires liberally)">
##INFO=<ID=DBSCAN_RB,Number=2,Type=Integer,Description="With --phasing both: DBSCAN allele lengths relative to reference, reported when DISCORDANT_LENGTH is set">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=RB,Number=2,Type=Integer,Description="Repeat length of the two alleles in bases relative to reference">
##FORMAT=<ID=FRB,Number=2,Type=Integer,Description="Full repeat length of the two alleles in bases">
##FORMAT=<ID=MRL,Number=2,Type=Integer,Description="Median read length of the two alleles' clusters in bases relative to reference, robust to a long length tail">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set identifier">
##FORMAT=<ID=SUP,Number=2,Type=Integer,Description="Read support per allele">
##FORMAT=<ID=SC,Number=2,Type=Integer,Description="Consensus score per allele">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  HG00271.hg38
chr1    1435798 .       TGGCGCGGAGCGGCGCGGAGCG  GCTGGCGCGGAGCGGCGCGGA,GCGGGCGCGCGCAGGA  .       .       END=1435818;STDEV=1,2     GT:RB:FRB:MRL:SUP:SC        1|2:1,-4:21,16:1,-4:18,6:63,41
chr1    57367044        .       AAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAATAAAT     AAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAATAAA,AAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAATAAA        .       .       END=57367125;STDEV=3,0 GT:RB:FRB:SUP:SC 1|2:-9,-34:72,47:17,12:216,141
```

## Tuning and hardcoded parameters

STRdust was developed while investigating the [pathogenic GOLGA8A repeat expansion](https://www.nature.com/articles/s41588-026-02537-7), and a [benchmark of repeat genotypers](https://www.biorxiv.org/content/10.64898/2026.02.28.708646v1) found it to be the most sensitive of the tested tools for detecting actual pathogenic expansions. Maximising that sensitivity is a priority, and several of the heuristics that support it are intentionally kept as **hardcoded constants** rather than command-line arguments. A parameter is only worth exposing if a user can tell how to tune it, and exposing all of these would inflate the option list without helping most users. Each is defined as a documented constant in the source so it can be changed (and promoted to a CLI argument) if real cohort experience shows a need:

| Parameter | Value | Where (source) | Tuning intuition |
|---|---|---|---|
| `EXPANSION_OUTLIER` min reads / ratio | ≥ 2 reads, > 2x longer allele | `src/genotype.rs` | lower either for more sensitivity to few-read / single-read expansions, at the cost of more artifact flags |
| `IMPRECISE_LENGTH` length CV | > 0.2 | `src/consensus.rs` | lower flags more loci as length-imprecise; raise to flag only the most spread-out |
| `DISCORDANT_LENGTH` ratio (`--phasing both`) | > 2x longer allele | `src/genotype.rs` | lower to surface smaller Ward/DBSCAN disagreements for QC review |
| DBSCAN neighbourhood radius (`eps`) | 0.4 | `src/phase_insertions.rs` | smaller = tighter/more clusters (more splitting, more reads dropped to noise → can undercall length-variable expansions); larger = looser clusters that can merge distinct alleles. Adjust in ±0.1 steps and watch `NCLUSTERS` |
| DBSCAN length weight | 0.3 | `src/phase_insertions.rs` | lower = composition dominates (keeps length-variable expansions together); higher = length matters more (separates same-motif alleles that differ only in length, behaving more like `ward`) |

If you have a specifically challenging repeat expansion where the defaults do not work well, we would be happy to work with you — please [open an issue](https://github.com/wdecoster/STRdust/issues) with details.

## Development

Contributor setup, the development workflow, testing (including the
network-dependent `TEST_PATHOGENIC_NETWORK` / `TEST_PATHOGENIC_FULL` tests), code
quality standards, and CI are documented in [CONTRIBUTING.md](CONTRIBUTING.md).

## CITATION

If you use this tool, please consider citing our [publication](https://genome.cshlp.org/content/early/2024/08/15/gr.279265.124).
