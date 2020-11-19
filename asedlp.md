---
title: Prediction of allele-specific expression uing RNA-seq results by deep learning
author: Zhenhua Zhang
date: Nov 19, 2020
output:
  html_document:
    toc: true
---

# TODO
- [x] Tensorboard to check the training
- [x] Multiple classes ROC curve and precision-recall curve
- [x] Check the feasibility of phasing BIOS genotypes, if possible, then I'll
  have more than 3K samples for training. It's possible to phase the BIOS data
  set.


# Introduction


# Methods (data analysis)

## Preprocessing and QC
`fastp` preprocessing

## Reads mapping
`STAR` two-pass mapping

## Remove reference bias
`WASP` pipeline to remove reference bias

## Phasing and haplotype

1. **How the input files look like?**  
**Genotypes**  
The genotypes are hosted on BOXY cluster.  
    - TriTyper: /groups/umcg-bios/prm02/projects/genotypes  # Imputed but not phased  
    - Plink: /groups/umcg-bios/prm02/projects/imputed  # Including raw and  
      imputed (unknown pipeline)  

2. **How to quality control?**  
    - `GenomicHarmonizer` to convert TriTyper format into vcf  
    - `shapeit4` phase the genome for each individual  

## Quantification of ASE
In-house scripts to construct regulating matrix and to estimate ASE effects  
    - Regulating sequence
    - ASE effects

## Deep learning model
`PyTorch` to implement a `CNN` model


# Results

## Basic statistics
1. Merge P-values by beta-binomial test for each gene.
2. Adjust P-values of multiple-test by FDR.
3. Assign genes with FDR < 0.05 as ASE genes for the given individuals.
4. Stacked bar plot to show the percentage of gene with ASE. Perhaps, multiple
   plots with different FDR threshold are more expressive.
5. A Manhattan plot to show P-values across the genome per chromosome.
6. Gene enrichment in different gene sets. E.g. GO annotations, pathways.

## Model training and evaluation

### Training
1. CNN model (Adam + Cross-Entropy-Loss + Drop-out)
2. Cross-validation (not yet), the drop-out is an more convenient way to
   diminish over-fit

### Evaluation
1. Accuracy
2. Precision
3. Recall
4. ROC-AUC curve

### Class activation map (CAM)
A way to show important pixels which determine the class for a given input
matrix. It combines the last features map in forward propagation and fully
connected layers that includes classification information. The integration are
represented by heatmap in which the color emphasizes the importance of pixels
with respect to classification.


# Expanding target cohort from GoNL (268) to BIOS (~3000)
In the pilot study to predict ASE, the results demonstrated that features in
upstream sequence of a gene are reliable predictors to discriminate ASE from
non-ASE genes. However, there are risks that the CNN model could be over-fitted
as the small sample size (i.e. 268 samples). To minimize the over-fitting
issue, more samples will be exploited in the model.


## Data

We firstly work on chr 17 ase it represents a chromosome with a minimal number
of parent-specific imprinted genes[^1].

## GRCh37 reference

## BIOS RNA-seq results

## GoNl genotypes and BIOS genotypes



## Reference
[^1]: The landscape of genomic imprinting across diverse adult human tissues.


