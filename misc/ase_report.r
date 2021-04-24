#!/usr/bin/env Rscript
library(dplyr)
library(ggplot2)

workdir <- "~/Documents/projects/wp_ase_dlp/outputs/aseReports"
setwd(workdir)

gene_list <- paste0(workdir, "/gene_ratio_gt_0.1_per_cohort.csv")
gene_dtfm <- read.csv(gene_list, header=TRUE)

cohort_comb <- gene_dtfm[, 'cohort'] %>%
    as.character() %>%
    unique() %>%
    combn(m=2) %>%
    t() %>%
    as.data.frame()

gene_dtfm %>%
    group_by(gene_id) %>%
    mutate(ncohort=length(cohort)) %>%
    ggplot() + geom_density(aes(x=ncohort))

