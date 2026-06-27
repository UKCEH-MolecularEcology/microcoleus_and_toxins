# Microcoleus Genome Abundance Workflow (Snakemake)

This workflow estimates the abundance of **Microcoleus genomes**
(defined by amino acid marker sequences in
`Dataset2.AA.microcoleus.and.GenBank.fa`) across assembled metagenomic
samples.

For each sample, the pipeline:

1.  Maps assembled contigs to Microcoleus AA markers using DIAMOND
    blastx
2.  Links contigs → Genome IDs
3.  Extracts read counts from existing BAM files
4.  Computes:
    -   Reads
    -   CPM (Counts Per Million)
    -   RPKM
    -   MeanDepth (length-weighted)
    -   Breadth (length-weighted coverage fraction)
5.  Produces:
    -   Per-sample abundance tables
    -   Genome × Sample abundance matrices

------------------------------------------------------------------------

## Input Requirements

### 1. Assembly directory

Defined in `config/config.yaml`:

``` yaml
assembly_dir: "assembly"
```

Directory must contain:

    A21_2022.fasta
    A21_2022.bam
    A27_2023.fasta
    A27_2023.bam
    ...

Naming requirement:

    <sample>.fasta  ↔  <sample>.bam

The sample name (basename without extension) must match exactly.

------------------------------------------------------------------------

### 2. Microcoleus marker FASTA

    Dataset2.AA.microcoleus.and.GenBank.fa

Defined in config:

``` yaml
aa_fasta: "Dataset2.AA.microcoleus.and.GenBank.fa"
```

Headers must contain genome IDs:

    >POL7_B2
    MAGGIYLIQDDDRL...

------------------------------------------------------------------------

## Software Requirements

The workflow uses tools from an existing conda environment:

    microcoleus

It must contain:

-   diamond
-   samtools
-   snakemake (optional)
-   awk
-   python

------------------------------------------------------------------------

## Running the Workflow

Activate environment containing Snakemake:

``` bash
conda activate snakemake
```

Run:

``` bash
snakemake   --configfile config/config.yaml   -s Snakefile   --cores 64   --jobs 8   -rp
```

Recommended flags:

-   `-r` → show reason
-   `-p` → print shell commands
-   `--show-failed-logs` → debugging
-   `-n` → dry-run

Example dry-run:

``` bash
snakemake -npr --cores 64
```

------------------------------------------------------------------------

## Parallelisation Strategy

Configured for HPC:

-   64 cores total
-   8 concurrent jobs

Per rule:

-   index_bam: 4 threads
-   contig2genome: 8 threads

------------------------------------------------------------------------

## Output Structure

    microcoleus_abundance_out/
    │
    ├── diamond_db/
    │   └── Dataset2_AA.dmnd
    │
    ├── tmp/
    │   ├── genome_ids.txt
    │   ├── bam_indexed/
    │   └── contig2genome/
    │
    ├── per_sample/
    │   ├── A21_2022.genome_abundance.tsv
    │   └── ...
    │
    └── matrices/
        ├── genome_Reads_matrix.tsv
        ├── genome_CPM_matrix.tsv
        ├── genome_RPKM_matrix.tsv
        ├── genome_MeanDepth_matrix.tsv
        └── genome_Breadth_matrix.tsv

------------------------------------------------------------------------

## Per-Sample Output Format

    GenomeID  Sample  Reads  CPM  RPKM  MeanDepth  Breadth  ContigsHit  SumContigLen_bp  TotalMappedReads

------------------------------------------------------------------------

## Matrices

Each matrix format:

    GenomeID    Sample1    Sample2    Sample3 ...

Available metrics:

-   genome_Reads_matrix.tsv
-   genome_CPM_matrix.tsv
-   genome_RPKM_matrix.tsv
-   genome_MeanDepth_matrix.tsv
-   genome_Breadth_matrix.tsv

------------------------------------------------------------------------

## Re-running

Force recomputation:

``` bash
snakemake --forceall --cores 64
```

Clean outputs:

``` bash
rm -rf microcoleus_abundance_out
```

------------------------------------------------------------------------

## Notes

-   This workflow uses `conda run -n microcoleus` inside rules.
-   Do NOT use `--use-conda` with this setup.
-   Designed for HPC (64 cores).
