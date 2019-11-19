# Allele-Specific Expression Deep Learning Predictor

## Mapping reads to reference genome
- Quality control of raw reads was done using fastp
- Clean reads were mapped to reference genome using Hisat2
- WASP pipeline was used to smooth the allelic mapping bias
- Indel realignment and base quality score calibration was done using GATK3

## Construct UTR sequence matrix

## Quantify ASE effects from RNA seq results

## Train deep learning(DL) model using TensorFlow

## Validate the DL model using GTEx


## Directory tree
/home/umcg-zzhang/Documents/projects/ASECausalSNPPrioritization
|-- buffers             Output files will be used many times or as breakpoints between input and final output.
|-- inputs              Any files not produced during current project.
|-- misc                Any files not grouped into other directories.
|-- outputs             Final output that will be published or shared.
|-- README.md           This file itself
|-- scripts             Any scripts use in current project
`-- temps               Temporary files


## FAQ
1. What does the Inaccessible in the FILTER field mean?
SNPs in "inaccessible' genome region, the inaccessible regions are enriched for novel SNVs; many are likely to be false positives.

2. How many individuals in GoNl cohort (only parents)
There are 499 individuals in total, 250 fathers, 248 mothers, 1 mothers (should be excluded in the anlysis) with haplotypes inferred from the children.

3. How many individuals(only parents) both in GoNl and BIOS
There are 397 individuals, including 134 fathers, 134 mathers and 129 children. **Therefore, we only have at most 268 samples to be analyzed.** Questions could be that: is the sample size large enough?

4. Which reference genome was used in the GoNl project? By which aligner?
Refrence genome: UCSC human reference genome build 37(`/home/cog/lfrancioli/resources/hg19/human_g1k_v37.fa`). 
Aligner: BWA(0.5.9-r16)

5. Which reference genome, gene annotation information, connonical indels, and snps will be used?
Reference genome: human_g1k_v37.fasta
Gene annotation info: ensembl GRCh37 release-75
Connonical INDELs:
Connonical SNPs:

6. How to eliminate the allelic mapping bias?
We have many options but we chose WASP to correct the allelic mapping bias.
Actually I had been hasitating between hard-mask [by Niek et.] and WASP pipepline. This [post](https://www.biostars.org/p/290455/) pushed me to choose the later. Refere to Lude's paper on an allelic expression analysis on GoNL and
BIOS cohort.

7. How to compile `snp2h5 / fasta2h5` on cluster?
> Compile snp2h5 is only needed to be done if you plan use snp2h5 or fast2h5
```
# 1. The following version of HDF5 functions in the compilation
$> module load HDF5/1.8.14-foss-2015b

# 2. Change `HDF_INSTALL=${HOME}/anoconda` to the path where HDF5 is installed

# 3. Optional
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/apps/software/

# 4. make
cd snp2h5
make
```

8. How to quantify allelic expression
- ASElux
    - No quality control
- STAR + WASP
    - With quality control but slow


9. Competitions
A paper using BIOS and GoNL cohort was published in this June (Lude is one of the authors).