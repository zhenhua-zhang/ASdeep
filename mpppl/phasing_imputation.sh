#!/bin/bash
set -Eeu -o pipefail

source /apps/modules/modules.bashrc
pjdir=~/Documents/projects/wp_ase_dlp

# Cohorts
# CODAM  NTR  PAN  RS

sub_cohort_id=CODAM
chr_id=chr22
thread=10
seed=31415

sbatch --time 0:29:0 \
    --mem 40G \
    --cpus-per-task ${thread} \
    --output %j-%u-beagle-${sub_cohort_id}-${chr_id}.log \
    --job-name beagle-${sub_cohort_id}-${chr_id} \
    <<EOF | cut -f4 -d' '
#!/bin/bash
set -Eeux -o pipefail

# Ensure Slurm settings is updated.
source /apps/modules/modules.bashrc

#
## Quality control by PLINK
#
module load plink
module list

plink --version

# Plink basic statistics, inlcuding heterogenity allele frequency,
# Hardy-Weinberg quibrilium individual and genotype missing rate.
# No --check-sex, since gender doesn't matter much.
plink --bfile ${sub_cohort_id} \
    --het \
    --freq \
    --hardy \
    --missing \
    --out ${sub_cohort_id}-${chr_id}-plink_stat

# Filter out low quality genotypes and encode the genotypes in VCF.
# --maf  Minor allele frequency
# --hwe  Hardy-Weinberg equilibrium
# --geno Missing call frequencies
plink --bfile ${sub_cohort_id} \
    --maf 5e-3 \
    --hwe 1e-3 \
    --geno 0.05 \
    --recode vcf bgz \
    --out ${sub_cohort_id}-${chr_id}


#
## Phase and imputation by Beagle
#

# Beagle version 5.1 which is powered by Java 1.8 or higher.
module load Java/1.8.0_45
module list
java -jar ~/tools/beagle/beagle.18May20.d20.jar | head -1

# Imputation and phasing.
java -jar ~/tools/beagle/beagle.18May20.d20.jar \
    gt=${pjdir}/buffers/${sub_cohort_id}-${chr_id}.vcf.gz \
    ref=${pjdir}/inputs/beagle/reference_panel/${chr_id}.1kg.phase3.v5a.vcf.gz \
    map=${pjdir}/inputs/beagle/genetic_map/plink.${chr_id}.GRCh37.map \
    seed=${seed} \
    nthreads=${thread} \
    out=${pjdir}/buffers/beagle-${sub_cohort_id}-${chr_id}
EOF
