# Allele-Specific Expression Deep Learning Predictor - asedlp

# Overview of the pipeline

![Overview-of-the-asedlp-pipeline.png]()

## Mapping reads to reference genome
- Quality control of raw reads was done using `fastp`
- Clean reads were mapped to reference genome using `STAR`
- `WASP` pipeline was used to smooth the allelic mapping bias
- `WASP` pipeline to count reads per haplotype

## Construct UTR sequence matrix

## Quantify ASE effects from RNA seq results

## Train deep learning model using PyTorch

## Validate the deep learning model using GTEx dataset

## Directory tree
Current dir: `~/Documents/projects/wp_ase_dlp`
```
.
├── README.md   This file itself
├── inputs      Any files not produced during current project.
├── buffers     Output files will be used many times or as breakpoints between input and final output.
├── outputs     Final output that will be published or shared.
├── scripts     Any scripts use in current project
├── misc        Any files not grouped into other directories.
└── temps       Temporary files
```


# Overview of dataset

## On BIOS VMs

## BIOS on UMCG clusters
There are four cohort with genotypes and RNA-seq results available. In total,
x,xxx individuals consists of the whole BIOS dataset.

### Path to genotypes (imputed by Michigan Imputation Server)
In total, there are 3,768 genotypes available. For each cohort, there are
| Cohort       | Samples | Path (Calculon) |
| :----------: | :-----: | :---:           |
| CODAM        | 561     | /groups/umcg-bios/prm02/projects/HRC_imputation/CODAM/results/unzipped
| NTR_GoNL     | 333     | /groups/umcg-bios/prm02/projects/HRC_imputation/NTR/GoNL/results/unzipped
| NTR_Aff6     | 1805    | /groups/umcg-bios/prm02/projects/HRC_imputation/NTR/Affy6/results/unzipped
| PAN          | 192     | /groups/umcg-bios/prm02/projects/HRC_imputation/PAN/results/unzipped
| RS           | 877     | /groups/umcg-bios/prm02/projects/HRC_imputation/RS/results/unzipped
| LLS_OminExpr | 236     |
| LLS_660Q     | 377     |
| LL           | 1134    |

**TODO**: Don't know the exact imputation settings yet.

### Path to RNA-seq files (Calculon / Gearshift?)
The RNA-seq FASTQ files could be found: `/groups/umcg-bios/prm02/rawdata/rnaseq`

### After matching genotype ids to RNA-seq ids
More info at: `~/Documents/projects/wp_ase_dlp/inputs/idmapping`
| Cohort_id    | Successful_count  | Schedualed_count  | Pass  |
| :----------: | :---------------: | :---------------: | :---: |
| CODAM        | 180               | 180               | Yes   |
| LLS_660Q     | 372               | 372               | Yes   |
| LLS_OminExpr | 236               | 236               | Yes   |
| NTR_Affy6    | 744               | 744               | Yes   |
| PAN          | 167               | 167               | Yes   |
| RS           | 698               | 698               | Yes   |
| gonl         | 273               | 273               | Yes   |
| LL           | 0                 | 407               | No    |
| Total        | 2670              | 3077              | No    |

**Some issues:**  
1. The LL deep are note imputed by Michigan imputation server yet.


## Geuvadis

### RNA-seq FASTQ files

### HD genotype chip data


## GTEx
We requested GTEx dataset (including 979 participants) by Project dbGap Accession Number: phs000424.v8.p2.c1

0. 

1. **RNA-seq data**
We downloaded BAM files for 922 whole blood samples (Whole Blood: Whole Blood) using Gen3-Client following the instruction at [GTEx v8 free egress instructions](https://anvilproject.org/learn/reference/gtex-v8-free-egress-instructions).
We used the newest (9th Dec. 2020) Linux version of Gen3-Client downloaded from [GitHub](https://github.com/uc-cdis/cdis-data-client/releases/download/2020.12/dataclient_linux.zip)
The dataset includes 618 males (median age xxx) and 304 females (median age xxx)

2. VCF files?


## FAQ
1. What does the Inaccessible in the FILTER field mean?
SNPs in "inaccessible' genome region, the inaccessible regions are enriched for
novel SNVs; many are likely to be false positives.

2. How many individuals in GoNl cohort (only parents)
There are 499 individuals in total, 250 fathers, 248 mothers, 1 mothers (should
be excluded in the anlysis) with haplotypes inferred from the children.

3. How many individuals(only parents) both in GoNl and BIOS
There are 397 individuals, including 134 fathers, 134 mathers and 129 children.  **Therefore, we
only have at most 268 samples to be analyzed.** Questions could be that: is the sample size large
enough?

4. Which reference genome was used in the GoNl project? By which aligner?
Refrence genome: UCSC human reference genome build 37 (/home/cog/lfrancioli/resources/hg19/human_g1k_v37.fa).
Aligner: BWA(0.5.9-r16)

5. Which reference genome, gene annotation information, connonical indels, and snps will be used?
Reference genome: human_g1k_v37.fasta
Gene annotation info: ensembl GRCh37 release-75
Connonical INDELs:
Connonical SNPs:

6. How to eliminate the allelic mapping bias?
We have many options but we chose WASP to correct the allelic mapping bias.
Actually I had been hasitating between hard-mask [by Niek et.] and WASP
pipepline. This [post](https://www.biostars.org/p/290455/) pushed me to choose
the later. Referee to Lude's paper on an allelic expression analysis on GoNL and
BIOS cohort.

7. How to compile `snp2h5 / fasta2h5` on cluster?
> Compile snp2h5 is only needed to be done if you plan use snp2h5 or fast2h5

```
# 1. The following version of HDF5 functions in the compilation
$> module load HDF5/1.8.14-foss-2015b

# 2. Change `HDF_INSTALL=${HOME}/anoconda` to the path where HDF5 is installed

# 3. Optional
$> export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/apps/software/

# 4. make
$> cd snp2h5
$> make
```

8. How to quantify allelic expression
- ASElux
    - No quality control
- STAR + WASP
    - With quality control but slow

9. Competitions
A paper using BIOS and GoNL cohort was published in this June (Lude is one of
the authors).

10. Remove duplicated reads?
In the WASP pipeline, the duplicated reads will be removed by rmdup_pe.py or
rmdup.py for paired-end or single-end reads. But for RNA-seq experiments, is it
reasonable to remove the duplicates?


11. Where the BIOS genotypes are?

```
# The raw genotypes in PLINK format were located at
/groups/umcg-bios/prm02/projects/imputed/CODAM/genotypes/b37
```



# The Outline of the paper

## Prediction of allele-specific expression using RNA-seq results by deep learning

## Abstract

## Introduction

## Methods (data analysis)

### Phasing and haplotype
#### How the input files look like
The genotypes are hosted on BOXY cluster.  
    - TriTyper: /groups/umcg-bios/prm02/projects/genotypes  # Imputed but not phased, these files are not convertable because of duplicated SNP ids
    - Plink: /groups/umcg-bios/prm02/projects/imputed  # Including raw and imputed (unknown pipeline)  

#### How to quality control
    - `GenomicHarmonizer` to convert TriTyper format into vcf
    - `shapeit4` phase the genome for each individual

### ASE quantification
#### Preprocessing and QC
`fastp` preprocessing

#### Reads mapping
`STAR` two-pass mapping

#### Remove reference bias
`WASP` pipeline to remove reference bias

#### Quantification of ASE
In-house scripts to construct regulating matrix and to estimate ASE effects  
    - Regulating sequence
    - ASE effects

Please note:
    1. MHC (HLA) region should be excluded.
    2. Current only SNPs should be included, because indels cause a different
    length of two allele.

### Deep learning model
`PyTorch` to implement a `CNN` model, using pre-constructed ResNext
architechture


## Results

### Basic statistics
1. Merge P-values by beta-binomial test for each gene.
2. Adjust P-values of multiple-test by FDR.
3. Assign genes with FDR < 0.05 as ASE genes for the given individuals.
4. Stacked bar plot to show the percentage of gene with ASE. Perhaps, multiple
   plots with different FDR threshold are more expressive.
5. A Manhattan plot to show P-values across the genome per chromosome.
6. Gene enrichment in different gene sets. E.g. GO annotations, pathways.

### Model training and evaluation

#### Training
1. CNN model (Adam + Cross-Entropy-Loss + Drop-out)
2. Cross-validation (not yet), the drop-out is an more convenient way to
   diminish over-fit

#### Evaluation
1. Accuracy
2. Precision
3. Recall
4. ROC-AUC curve

#### Class activation map (CAM)
A way to show important pixels which determine the class for a given input
matrix. It combines the last features map in forward propagation and fully
connected layers that includes classification information. The integration are
represented by heatmap in which the color emphasizes the importance of pixels
with respect to classification.

### Expanding target cohort from GoNL (268) to BIOS (~3000)
In the pilot study to predict ASE, the results demonstrated that features in
upstream sequence of a gene are reliable predictors to discriminate ASE from
non-ASE genes. However, there are risks that the CNN model could be over-fitted
as the small sample size (i.e. 268 samples). To minimize the over-fitting
issue, more samples will be exploited in the model.

### Data

We firstly work on chr 17 ase it represents a chromosome with a minimal number
of parent-specific imprinted genes[^1].

### GRCh37 reference

### BIOS RNA-seq results

### GoNl genotypes and BIOS genotypes

### Pre-process

#### Match RNA-seq ID to DNA array ID

## Reference
[^1]: The landscape of genomic imprinting across diverse adult human tissues.

<!-- vim: set nospell: -->
