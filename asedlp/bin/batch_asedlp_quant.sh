#!/bin/bash

# Author     : Zhenhua Zhang
# Email      : zhenhua.zhang217@gmail.com
# License    : MIT
# Create date: Mon 09 Mar 2020 09:22:50 AM CET
# Last update: Mon 30 Mar 2020 04:31:20 PM CEST

# One time script
pjdir=/groups/umcg-bios/tmp04/umcg-zzhang/projects/ASEDLP

while read fastq_id; do
    sbatch \
        --time 0:29:0 \
        --cpus 1 \
        --mem 5G \
        --array 1-22 \
        --job-name asedlp_quant \
        --output $pjdir/workdir/logdir/%A_%a-%u-${fastq_id}_asedlp_quant.log \
        asedlp_quant.sh \
        -w $pjdir/workdir \
        -i $fastq_id \
        -a $pjdir/inputs/Ensembl_references/Homo_sapiens.GRCh37.75.gtf \
        -s $pjdir/misc/freeze2_GoNL_related_GTE_30092016_QCpassed.csv \
        -G $pjdir/inputs/Ensembl_references \
        -v $pjdir/scripts/.env/bin/activate
done < 
