#!/bin/bash
#SBATCH --mem=5G
#SBATCH --time=00:25:00
#SBATCH --cpus-per-task=5
#SBATCH --output=%A_%a_geuvadis_genotype_qc.log
#SBATCH --job-name=geuvadis_genotype_qc

chrid=${SLURM_ARRAY_TASK_ID:=22}

pjdir=~/Documents/projects/wp_ase_dlp
smpfile=$pjdir/inputs/geuvadis/samples/igsr_Geuvadis-mRNA_Phase3_sample-name.txt
test ! -e $smpfile && echo "Failed to find file: ${LINENO}" && exit -1

module purge
module load BCFtools
module list

wkdir=$pjdir/inputs/geuvadis/genotypes
for chr in {1..22}; do
    ipvcf=$wkdir/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.deploid.snp.vcf.gz
    opvcf=${ipvcf/ALL/Geuvadis}

    bcftools view \
        --threads 5 \
        --output-type z \
        --min-alleles 2 \
        --max-alleles 2 \
        --force-samples \
        --output-file $opvcf \
        --samples-file $smpfile \
        --include 'TYPE=="snp" && INFO/EUR_AF >= 0.0005' \
        $ipvcf

    bcftools index --threads 5 --force $opvcf
done

# Using 1..22 to keep variants sorted by chromosome
bcftools concat \
    --allow-overlaps \
    --output $wkdir/Geuvadis.all.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.deploid.snp.vcf.gz \
    --output-type z \
    --threads 5 \
    $wkdir/Geuvadis.chr{1..22}*.vcf.gz
rm -fr $wkdir/Geuvadis.chr{1..22}*.vcf.gz

bcftools index \
    -t \
    $wkdir/Geuvadis.all.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.deploid.snp.vcf.gz


echo "Done"
