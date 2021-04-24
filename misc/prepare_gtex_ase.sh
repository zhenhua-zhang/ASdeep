#!/bin/bash
projdir=~/Documents/projects/wp_ase_dlp
iptdir=$projdir/inputs
optdir=$projdir/inputs

ase_count_path=$iptdir/GTEx/ase/whlbld
liftover_matrix=$iptdir/GTEx/misc/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz

for ase_table in $ase_count_path/*.ase_table.tsv.gz; do
    echo ${ase_table}
    csvtk join -d$'\t' -f'VARIANT_ID;variant_id' $ase_table $liftover_matrix \
        | grep -ie TISSUE_ID -e WHLBLD \
        | csvtk cut -I -f31,30,9,10,11 \
        | csvtk sep -I -f1 -s'_' -n'config,position,refAllele,altAllele' --drop \
        | csvtk cut -I -f6,7,2,8,9,3,4,5 \
        | csvtk rename -I -f3,6,7,8 -n variantID,refCount,altCount,totalCount \
        > ${ase_table/.tsv.gz/.whlbld.csv}
done
