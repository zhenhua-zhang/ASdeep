#!/bin/bash

chromid=${1}
fastqid=${2}
sampleid=${3}

chromid=${chromid:=1}
fastqid=${fastqid:=AC47H5ACXX-3-16}
sampleid=${sampleid:=gonl-101b}

source .env/bin/activate
./bin/asedlp quant \
    --haplotypes ../../workdir/snph5db/${chromid}/haplotype.h5 \
    --snp-tab ../../workdir/snph5db/${chromid}/snps_tab.h5 \
    --snp-index ../../workdir/snph5db/${chromid}/snps_index.h5 \
    --sequence-tab ../../workdir/fastah5db/chr${chromid}.h5 \
    --ref-read-counts ../../workdir/optdir/${fastqid}/waspOptdir/perChrom/${chromid}/${fastqid}_${chromid}.refAlleleCounts.h5 \
    --alt-read-counts ../../workdir/optdir/${fastqid}/waspOptdir/perChrom/${chromid}/${fastqid}_${chromid}.altAlleleCounts.h5 \
    --genome-annot ../../inputs/Ensembl_references/Homo_sapiens.GRCh37.75.gtf \
    --sample-id $sampleid \
    --gene-ids  ENSG00000008130 ENSG00000069424 ENSG00000142634 ENSG00000162511 ENSG00000134686 ENSG00000162434 ENSG00000143390 ENSG00000160710 ENSG00000158869 ENSG00000162747 ENSG00000188404 ENSG00000235750 ENSG00000116701
