# Introduction


# Pipeline

## `fastp`
- Quality control, trimming etc.

## Mask reference genome
- SNPs with MAF > 0.01??


## `STAR`
- Build index (once)
- Mapping with two-pass mapping strategy.
- Genome `human_g1k_v37`, annotations `Homo_sapiens.GRCh37.75`
STAR is slow, we need a trade-off between quality and speed. 

## `WASP`
WASP use Python, therefore multiple cores are meaningless. To save my priority, separated script
is recommended.

A strategy here is:
	1. Two-pass strategy mapping FASTQ files to the reference genome.
	2. Split mapped reads by chromosome.
	3. Apply WASP pipeline by chromosome, meanwhile using shared memory of genome index?
	4. The job dependency can be controlled by `--dependency` of `sbatch`. Refer to [sbatch dependency](https://hpc.nih.gov/docs/job_dependencies.html)

- SNP to HDF5 (once)
    check if there is an HDF5 database then create an environmental variable, otherwise construct the database.
- remap

## `phASER`
[Step-by-step](https://stephanecastel.wordpress.com/2017/02/15/how-to-generate-ase-data-with-phaser/)
Remove HLA gene as there's a lot repeats in the gene.

## Other ideas
Quantification without mapping to the reference genome and transcriptome annotations
`Salmon`, `Kallisto`, and `Sailfish` are tools to quantify gene expression without mapping reads to the reference genome.

## Other options
1. Constructing parental genome by `AlleleSeq`, then mapping RNA-seq result to parental genomes by junction aligner e.g. STAR
2. Using `GSNAP` which has a variant-aware alignment mode.
3. Another tool is `AlleleCounter` which belongs to the same author of `phASER`

