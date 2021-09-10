#!/usr/bin/env Rscript
#
## A utility script to check the enrichment of identified ASE genes
#

library(dplyr)
library(ggplot2)
library(enrichplot)
library(org.Hs.eg.db)
library(clusterProfiler)

rm(list=ls())

min_ratio <- 0.1
cohort <- "Geuvadis"

wkdir <- paste("~/Documents/projects/wp_ase_dlp/outputs/aseReports", cohort, sep="/")
setwd(wkdir)

#
## Basic statistical analysis
#


#
## Estimate
#

file_path <- paste(wkdir, "asedlp_report_exMHC-pvm_sorted_filtered_adj.csv", sep="/")
pval_mtx <- read.csv(file_path, header=TRUE) %>%
    filter(ratio >= min_ratio) %>%
    as.data.frame()

cat("In total, there are ", nrow(pval_mtx),
    " genes with at least ", min_ratio * 100,
    "% individuals show allelic expression.\n", sep="")


#
## GO terms
#
ont_type="BP"
ergo <- enrichGO(gene=as.character(pval_mtx$X),
                 OrgDb=org.Hs.eg.db,
                 keyType="ENSEMBL",
                 ont=ont_type,
                 pAdjustMethod="BH",
                 pvalueCutoff=0.01,
                 qvalueCutoff=0.05,
                 readable=TRUE)

show_cat <- 20
h <- show_cat / 4 + 1
w <- max(nchar(head(ergo@result$Description, show_cat))) / 15 + show_cat / 3 + 1
ggsave(paste(ont_type, min_ratio, "dotplot.pdf", sep="-"),
       plot=enrichplot::dotplot(ergo, showCategory=show_cat),
       width=w,
       height=h)

ggsave(paste(ont_type, min_ratio, "emapplot.pdf", sep="-"),
       plot=enrichplot::emapplot(ergo, showCategory=show_cat),
       width=10,
       height=10)

ggsave(paste(ont_type, min_ratio, "upsetplot.pdf", sep="-"),
       plot=enrichplot::upsetplot(ergo, n=show_cat),
       width=w,
       height=h * 2)


# 
## KEGG pathways
#
if (FALSE) {
    ont_type <- 'KEGGpathway'
    erkegg <- enrichKEGG(gene=pval_mtx$X,
                         organism = "hsa", 
                         keyType="ENSEMBL",
                         pvalueCutoff = 0.05, 
                         pAdjustMethod = "BH", 
                         minGSSize = 10, 
                         maxGSSize = 500, 
                         qvalueCutoff = 0.2, 
                         use_internal_data = FALSE)

    h <- show_cat / 4 + 1
    w <- max(nchar(head(erkegg@result$Description, show_cat))) / 15 + show_cat / 3 + 1
    ggsave(paste(ont_type, min_ratio, "dotplot.pdf", sep="-"),
           plot=dotplot(erkegg, showCategory=show_cat),
           width=w,
           height=h)

    ggsave(paste(ont_type, min_ratio, "emapplot.pdf", sep="-"),
           plot=emapplot(erkegg, showCategory=show_cat),
           width=10,
           height=10)

    ggsave(paste(ont_type, min_ratio, "upsetplot.pdf", sep="-"),
           plot=upsetplot(erkegg, n=show_cat),
           width=w,
           height=h * 2)

}
