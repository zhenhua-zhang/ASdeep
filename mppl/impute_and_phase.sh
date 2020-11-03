#!/bin/bash
#
# File name : impute_and_phasing.sh
# Author    : Zhenhua Zhang
# E-mail    : zhenhua.zhang217@gmail.com
# Created   : Tue 30 Jun 2020 04:49:32 PM CEST
# Version   : V 0.1.0
# License   : MIT
#

# Here to find Plink format genotypes: /groups/umcg-bios/prm02/projects/imputed/
# The genotype has been imputed and phased(?), but it's still possible to run
# the pipeline from raw data. If so, the ID between samples should be rematched.

set -Eeux -o pipefail

pjdir=${HOME}/Documents/projects/wp_ase_dlp
ipdir=${pjdir}/inputs
opdir=${pjdir}/outputs
bfdir=${pjdir}/buffers

memory=16G
threads=5
chr=chr1

# QC
plink_check_err_jobid=$(sbatch --time 0:29:0 \
    --mem ${memory} \
    --cpus-per-task ${threads} \
    --output %j-%u-plink_check_err-${chr}.log \
    --job-name plink_check_err-${chr} <<EOF | cut -f4 -d$' '
#!/bin/bash
set -Eeux -o pipefail

source /apps/modules/modules.bashrc
pjdir=
EOF
)


# Phasing
shapeit_jobid=$(sbatch --time 1:59:0 \
    --mem ${memory} \
    --cpus-per-task ${threads} \
    --output %j-%u-shapeit-${chr}.log \
    --dependency afterok:"${plink_check_err_jobid}" \
    --job-name shapeit-${chr} <<EOF | cut -f4 -d$' '
#!/bin/bash
set -Eeux -o pipefail

source /apps/modules/modules.bashrc
module load shapeit
shapeit \
    -B ${chr} \
    -M ${ipdir}/1kg_p3/gmap/${chr}.b37.gmap.gz \
    -T ${threads} \
    -O ${bfdir}/${chr}-shapeit \
    --seed 31415
EOF
)


# Impute
impute_jobid=$(sbatch --time 0:20:0 \
    --mem 16G \
    --cpus-per-task 8 \
    --output %j-%u-impute2-${chr}.log \
    --dependency=afterok:"${shapeit_jobid}" \
    --job-name impute2-${chr} <<EOF | cut -f4 -d$' '
#!/bin/bash
set -Eeux -o pipefail

source /apps/modules/modules.bashrc
module load IMPUTE2

impute2 \
    -use_prephased_g \
    -m ${ipdir}/1kg_p3/1000GP_Phase3/genetic_map_${chr}.txt \
    -h ${ipdir}/1kg_p3/1000GP_Phase3/1000GP_Phase3_${chr}.hap.gz \
    -l ${ipdir}/1kg_p3/1000GP_Phase3/1000GP_Phase3_${chr}.legend.gz \
    -known_haps_g ${pjdir}/buffers/${chr}-phased.haps \
    -Ne 20000 \
    -int 1 5e6 \
    -o ${bfdir}/${chr}-shapit-impute2 \
    -phase

    # -strand_g Strand alignments is not usually necessary when you just want to
    # phase a dataset, but it is important when that dataset will be combined
    # with a reference panel in a down stream analysis
EOF
)


# Convert to VCF format.
sbatch --time 0:20:0 \
    --mem 4G \
    --cpus-per-task 8 \
    --output %j-%u-bcftools_convert-${chr}.log \
    --dependency afterok:"${impute_jobid}" \
    --job-name bcftools_convert-${chr} <<EOF
#!/bin/bash
set -Eeux -o pipefail

source /apps/modules/modules.bashrc
module load BCFtools

bcftools convert \
    -G ${bfdir}/${chr}-shapit-impute2
EOF
