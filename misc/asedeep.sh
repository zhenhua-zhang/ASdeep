#!/bin/bash
wkdir=/home/umcg-zzhang/Documents/projects/wp_ase_dlp
source $wkdir/scripts/.env/bin/activate

asdeep makedb \
  -g $wkdir/../wp_reference/Ensembl/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
  -v $wkdir/temps/dpp9.vcf.gz \
  -i $wkdir/temps/dpp9.bed.gz \
  -m $wkdir/temps/dpp9.csv \
  -n 25000 \
  -o $wkdir/temps/dpp9

asdeep train -d $wkdir/temps/dpp9/DPP9.h5

asdeep predict \
  -t DC IG GS \
  -d $wkdir/temps/dpp9/DPP9.h5 \
  -s NA20818 \
  -m $wkdir/temps/asedeep_model.pth \
  -f png
