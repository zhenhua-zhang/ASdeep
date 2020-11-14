#!/bin/bash
#
# File name : phasing_imputation.sh
# Author    : Zhenhua Zhang
# E-mail    : zhenhua.zhang217@gmail.com
# Created   : Tue 30 Jun 2020 04:49:32 PM CEST
# Version   : V 0.1.0
# License   : MIT
#

set -Eeu -o pipefail

source /apps/modules/modules.bashrc
pjdir=/groups/umcg-bios/tmp03/users/umcg-zzhang/projects/wp_ase_dlp

# Sub dirs
plink_dir=$pjdir/outputs/phasingAndImputation/plink
beagle_dir=$pjdir/outputs/phasingAndImputation/beagle
log_dir=$pjdir/outputs/phasingAndImputation/logs
mkdir -p $plink_dir $beagle_dir $log_dir

# Cohorts: CODAM  NTR  PAN  RS
# Sample information
cohort_id=CODAM
chr_id=chr22

# Beagle settings
thread=5
seed=31415

# Sbatch settings
time=0:29:0
mem=40G
log=$log_dir/%j-%u-phasingAndImputation-$cohort_id-$chr_id.log
job=phasingAndImputation-$cohort_id-$chr_id

sbatch --time $time \
    --mem $mem \
    --output $log \
    --job-name $job \
    --cpus-per-task ${thread} \
    <<EOF | cut -f4 -d' '
#!/bin/bash
set -Eeu -o pipefail

# Ensure Slurm settings is updated.
source /apps/modules/modules.bashrc

#
## Quality control by PLINK
#
echo --------------------------------------------------------------------------
module purge
module load plink
module list

# PLINK version: PLINK v1.90b3.32 64-bit (24 Feb 2016)
echo --------------------------------------------------------------------------
plink --version

# Plink basic statistics, inlcuding heterogenity, allele frequency,
# Hardy-Weinberg quibrilium individual and genotype missing rate.
# No --check-sex, since gender doesn't matter much.
plink_bfile=$pjdir/buffers/$chr_id
plink_vfile=$plink_dir/plink-$cohort_id-$chr_id
echo --------------------------------------------------------------------------
plink --bfile \$plink_bfile \
    --het \
    --freq \
    --hardy \
    --missing \
    --out $plink_dir/plink-$cohort_id-$chr_id-plink_stat

# Filter out low quality genotypes and encode the genotypes in VCF.
# --maf  Minor allele frequency
# --hwe  Hardy-Weinberg equilibrium
# --geno Missing call frequencies
echo --------------------------------------------------------------------------
plink --bfile \$plink_bfile \
    --maf 5e-3 \
    --hwe 1e-4 \
    --geno 0.05 \
    --recode vcf bgz \
    --out \$plink_vfile


#
## Phasing and imputation by Beagle
#
# Beagle version 5.1 which is powered by Java 1.8 or higher.
echo --------------------------------------------------------------------------
module purge
module load Java/1.8.0_45
module list

# Beagle version: beagle.18May20.d20.jar (version 5.1)
echo --------------------------------------------------------------------------
java -jar ~/tools/beagle/beagle.18May20.d20.jar | head -1

# Imputation and phasing.
beagle_vgt=\$plink_vfile.vcf.gz
beagle_ref=$pjdir/inputs/beagle/reference_panel/$chr_id.1kg.phase3.v5a.vcf.gz
beagle_map=$pjdir/inputs/beagle/genetic_maps/plink.$chr_id.GRCh37.map
beagle_out=$beagle_dir/beagle-$cohort_id-$chr_id

echo --------------------------------------------------------------------------
java -jar ~/tools/beagle/beagle.18May20.d20.jar \
    gt=\$beagle_vgt \
    ref=\$beagle_ref \
    map=\$beagle_map \
    seed=$seed \
    nthreads=$thread \
    out=\$beagle_out
EOF
