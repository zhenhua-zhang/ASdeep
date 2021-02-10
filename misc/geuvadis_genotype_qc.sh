#!/bin/bash
#SBATCH --mem=5G
#SBATCH --array=1-22
#SBATCH --time=00:25:00
#SBATCH --cpus-per-task=5
#SBATCH --output=%A_%a_geuvadis_genotype_qc.log
#SBATCH --job-name=geuvadis_genotype_qc

chrid=${SLURM_ARRAY_TASK_ID:=22}

pjdir=~/Documents/projects/wp_ase_dlp
smpfile=$pjdir/inputs/geuvadis/samples/igsr_Geuvadis-mRNA_Phase3_sample-name.txt
test ! -e $smpfile && echo "Failed to find file: ${LINENO}" && exit -1

ipvcf=$pjdir/inputs/geuvadis/genotypes/ALL.chr$chrid.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
test ! -e $ipvcf && echo "Failed to find file: ${LINENO}" && exit -1

opvcf=${ipvcf/ALL/Geuvadis}

module purge
module load BCFtools
module list

bcftools view \
    --threads 2 \
    --output-type z \
    --min-alleles 2 \
    --max-alleles 2 \
    --force-samples \
    --output-file $opvcf \
    --samples-file $smpfile \
    --include 'TYPE=="snp" && INFO/EUR_AF >= 0.0005' \
    $ipvcf

bcftools index --threads 2 --force $opvcf

echo "Done"
