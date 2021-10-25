#!/usr/bin/env Rscript
library(data.table)
library(tidyverse)
library(liftOver)
library(magrittr)
library(biomaRt)

wkdir <- "../../inputs/Ensembl_references"

# Some code to explore the database
# listDatasets(ensembl) %>% head()
# searchAttributes(mart=ensembl, pattern="feature_page")
# listFilters(mart=ensembl) %>% fwrite("../../temps/biomart_filters.txt")
# listAttributes(mart=ensembl) %>% fwrite("../../temps/biomart_attrs.txt")

# We first fetch canonical and protein coding transcripts on autosomes from the
# GRCh38 build.
ensembl_38 <- useEnsembl(biomart="genes", dataset="hsapiens_gene_ensembl")

# NOTE: The `filters` should have the same order with `values`
flt_lt <- c("transcript_is_canonical", "biotype", "chromosome_name")
val_lt <- list(transcript_is_canonical=T, biotype="protein_coding",
               chromosome_name=as.character(1:22))
geneset_38 <- getBM(attributes="ensembl_transcript_id", filters=flt_lt,
                    values=val_lt, mart=ensembl_38)

# Then, we fetch the information of transcripts with the obtained GRCh38
# Ensembl IDs from GRCh37 build.
transcript_id_lst <- geneset_38$ensembl_transcript_id
ensembl_37 <- useEnsembl(biomart="genes", dataset="hsapiens_gene_ensembl",
                         version="GRCh37")

attrs <- c("chromosome_name", "exon_chrom_start", "exon_chrom_end",
           "ensembl_gene_id", "ensembl_transcript_id", "ensembl_exon_id",
           "external_gene_name", "external_gene_source", "strand",
           "gene_biotype")

geneset_37 <- getBM(attributes=attrs,
                    filters=c("ensembl_transcript_id", "chromosome_name"),
                    values=list(ensembl_transcript_id=transcript_id_lst,
                                chromosome_name=as.character(1:22)),
                    mart=ensembl_37)

geneset_37_gtf <- geneset_37 %>%
  dplyr::mutate(attr=paste0("gene_id \"", ensembl_gene_id, "\"; ",
                            "transcript_id \"", ensembl_transcript_id, "\"; ",
                            "exon_id \"", ensembl_exon_id, "\"; ",
                            "gene_name \"", external_gene_name, "\"; ",
                            "gene_source \"", external_gene_source, "\"; ",
                            "biotype \"", gene_biotype, "\";"),
                feature="exon",
                score=".",
                frame=".",
                strand=if_else(strand=="1", "+", "-")) %>%
  dplyr::rename("seqname"="chromosome_name", "start"="exon_chrom_start",
                "end"="exon_chrom_end", "source"="external_gene_source") %>%
  dplyr::select(one_of(c("seqname", "source", "feature", "start",
                         "end", "score", "strand", "frame", "attr"))) %>%
  dplyr::arrange(seqname, start, end)

opath <- str_glue("{wkdir}/Homo_sapiens.GRCh37.protein_coding.canonical.gtf")
geneset_37_gtf %>% fwrite(opath, sep="\t", col.names=F, row.names=F, quote=F)

