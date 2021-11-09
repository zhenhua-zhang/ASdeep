#!/bin/bash
set -Eeu -o pipefail
#
## Utilities {
#
# Error information, exit -1
err() {
  echo -e "[E]: $1" >&2 && exit -1
}

# Warning information
wrn() {
  echo -e "[W]: $*" >&2
}

# General information
msg() {
  echo -e "[I]: $*"
}
# } // Utilities


wkdir=~/Documents/projects/wp_ase_dlp
idmap_file=$wkdir/inputs/idmapping/BIOS-genotype-rnaseq-ids_with-LL_usable_20210105.txt
gff_file=$wkdir/inputs/Gencode/gencode.annotation.ensembl_canonical.protein_coding.GRCh37.gff3.gz

# Prepare read counts
for rcfile in $(ls $wkdir/temps/test_readcounts); do
  rc_name=$(basename $wkdir/temps/test_readcounts/$rcfile)
  msg "Converting $rc_name"

  rna_id=$(cut -f1 -d_ <<<$rc_name)
  var_id=$(grep -w $rna_id $idmap_file | cut -f2)
  awk -F, 'NR>1 {
      meta=$3;
      for(i=4;i<NF;i++) {meta=meta";"$i};
      meta=meta";"$NF;
      print $1"\t"$2-1"\t"$2"\t"meta}' $wkdir/temps/test_readcounts/$rcfile \
    | sort -k1,1g -k2,2g -k3,3g \
    | bgzip > $wkdir/temps/test_readcounts/$var_id.bed.gz

  tabix $wkdir/temps/test_readcounts/$var_id.bed.gz
done

# Prepare genotypes
# PLBD1: 12:14656595-14720817
# PYGB: 20:25228721-25278648
# Up/down stream 25000 bp
if [[ $(hostname) =~ gearshift ]]; then
  PLDB1=12:14631595-14745818
  PYGB=20:25203721-25303648

  bcftools view $wkdir/outputs/phasing/all-PAN-singleAlt.vcf.gz $PLDB1 $PYGB \
    | bgzip > $wkdir/temps/test.vcf.gz

  tabix $wkdir/temps/test.vcf.gz
else
  err "Genotypes are only available on the cluster (Gearshift)!"
fi

# Prepare GTF file
grep -w -e PLDB1 -e PYGB $gff_file \
  | bgzip > $wkdir/temps/test.gtf.gz
tabix $wkdir/temps/test.gtf.gz

# Test all sub-command
source $wkdir/scripts/.env/bin/activate

# Infer allelic imbalance
for rc_file in $wkdir/temps/test_readcounts/*.bed.gz; do
  sample_id=$(cut -f1 -d. <<<$(basename $rc_file))
  asdeep inferai \
    -v $wkdir/temps/test.vcf.gz \
    -i $wkdir/temps/test.gff.gz \
    -c $rc_file \
    -s $sample_id \
    -o $wkdir/temps/test_ase/$sample_id.ASE.csv
done

for per_ase_report in $wkdir/temps/test_ase/*.csv; do
  sample_id=$(basename $per_ase_report | cut -f1 -d.)

  while IFS="," read -a line; do
    gene_id=${line[0]}
    if [[ $gene_id =~ "gene_id" ]]; then continue; fi

    p_val=${line[2]}
    if [[ -z $p_val ]]; then continue; fi

    p_val_lf=${line[4]}
    p_val_rt=${line[5]}
    prone_a1=$(echo $p_val_rt'<0.5' | bc -l)
    prone_a2=$(echo $p_val_lf'>0.5' | bc -l)

    if [[ $prone_a1 -eq 1 ]]; then
      direction=-1
    elif [[ $prone_a2 -eq 1 ]]; then
      direction=1
    else
      direction=0
    fi

    opt_file=$wkdir/temps/$gene_id".ASE.csv"
    if [[ ! -e $opt_file ]]; then
      echo "sample_id,ASE" > $opt_file
    fi

    echo "$sample_id,$direction" >> $opt_file
  done < $per_ase_report
done

asdeep makedb \
  -g $wkdir/inputs/Ensembl_references/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
  -v $wkdir/temps/test.vcf.gz \
  -i $wkdir/temps/test.bed.gz \
  -m $wkdir/temps/test.csv \
  -n 32700 \
  -o $wkdir/temps/test_hbcdb

asdeep train \
  -d $wkdir/temps/test_hbcdb/PLBD1.h5 \
  -m $wkdir/temps/test_model/PLBD1.asdeep_model.pth \
  -L $wkdir/temps/test_model/PLBD1.asdeep_logs

asdeep train \
  -d $wkdir/temps/test_hbcdb/PYGB.h5 \
  -m $wkdir/temps/test_model/PYGB.asdeep_model.pth \
  -L $wkdir/temps/test_model/PYGB.asdeep_logs

asdeep predict \
  -t IG GS \
  -s ALS05219_ALS05219 \
  -d $wkdir/temps/test_hbcdb/PLBD1.h5 \
  -m $wkdir/temps/test_model/PLBD1.asdeep_model.pth \
  -f png

