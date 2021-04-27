#!/bin/bash
#
# File name : prepare_bios_vars.sh
# Author    : Zhenhua Zhang
# E-mail    : zhenhua.zhang217@gmail.com
# Created   : Tue 30 Jun 2020 04:49:32 PM CEST
# Version   : V 0.1.0
# License   : MIT
#

if [[ -e /apps/modules/modules.bashrc ]]; then
    source /apps/modules/modules.bashrc
fi

set -Eeu -o pipefail

# Cohorts: CODAM  NTR  PAN  RS LL
# Sample information
cohort="${1:?Error: You should give the cohort id!}"

# Sub dirs
#pjdir=/Documents/projects/wp_ase_dlp
pjdir=/groups/umcg-bios/tmp01/users/umcg-zzhang/projects/wp_ase_dlp
if [[ ! -d $pjdir/inputs/BIOS_genotypes/$cohort ]]; then
    echo "Failed to find $pjdir/inputs/BIOS_genotypes/$cohort! Exit"
    exit
fi

bcftools_dir=$pjdir/outputs/phasing_v2/$cohort/bcftools
conformGt_dir=$pjdir/outputs/phasing_v2/$cohort/conformGt
beagle_dir=$pjdir/outputs/phasing_v2/$cohort/beagle
log_dir=$pjdir/outputs/phasing_v2/$cohort/logs

mkdir -p $bcftools_dir $conformGt_dir $beagle_dir $log_dir

debug=true
if [[ $debug == true ]]; then
    runCpus=2
    runMems=15G
    runTime=0:9:0
    runArray=1-1
    set_flags=-Eeux
else
    runCpus=2
    runMems=15G
    runTime=19:59:0
    runArray=1-22
    set_flags=-Eeu
fi

# Beagle settings
seed=31415


sbatch \
    --mem=$runMems \
    --time=$runTime \
    --array=$runArray \
    --cpus-per-task=$runCpus \
    --job-name=phasing-$cohort \
    --output=$log_dir/%A-%a_%u_phasing_$cohort.log \
    <<EOF | cut -f4 -d' '
#!/bin/bash

# Ensure Slurm settings is updated.
if [[ -e /apps/modules/modules.bashrc ]]; then
    source /apps/modules/modules.bashrc
fi

set $set_flags -o pipefail

chr_id=\${SLURM_ARRAY_TASK_ID:=1}

module purge && module load BCFtools/1.9-foss-2018b && module list

# Input: upstream
# The imputed files are located /groups/umcg-bios/prm02/projects/HRC_imputation,
# which are imputed using Michigan Imputation Server.
imputed_gntp=$pjdir/inputs/BIOS_genotypes/$cohort/chr\${chr_id}.dose.vcf.gz

# Output: downstream
imputed_filtered_gntp=$bcftools_dir/\$(basename \${imputed_gntp/vcf.gz/filtered.vcf.gz})

# R2 refer to Xiaojing's QC on 200HIV genotypes (>=0.3);
# MAF is much smaller than the one used in GWAS analysis (>=0.1)
bcftools view \
    --include 'TYPE=="snp" && INFO/MAF >= 0.0005 && INFO/R2 >= 0.3' \
    --threads $runCpus \
    --output-type z \
    --output-file \$imputed_filtered_gntp \
    \$imputed_gntp


# Beagle version 5.1 which is powered by Java 1.8 or higher.
module purge && module load Java/1.8.0_144 && module list

# Input: reference
ref_panel_gntp=$pjdir/inputs/beagle/reference_panel/chr\${chr_id}.1kg.phase3.v5a.vcf.gz

# Input: upstream
# imputed_filtered_gntp

# Output: downstream
imputed_filtered_adj_gntp=$conformGt_dir/\$(basename \${imputed_filtered_gntp/vcf.gz/adjByRef})

java -jar ~/tools/beagle/conform-gt.24May16.cee.jar \
    gt=\$imputed_filtered_gntp \
    ref=\$ref_panel_gntp \
    chrom=\$chr_id \
    match=POS \
    out=\$imputed_filtered_adj_gntp

# Phasing.
# Input: reference map and genotype.
ref_panel_map=$pjdir/inputs/beagle/genetic_maps/plink.chr\${chr_id}.GRCh37.map

# Output: final
imputed_filtered_adj_phased_gntp=$beagle_dir/beagle-$cohort-chr\${chr_id}

# TODO: A blacklist of samples without RNA-seq results.

java -jar ~/tools/beagle/beagle.18May20.d20.jar \
    gt=\$imputed_filtered_adj_gntp".vcf.gz" \
    ref=\$ref_panel_gntp \
    map=\$ref_panel_map \
    impute=false \
    seed=$seed \
    nthreads=$runCpus \
    out=\$imputed_filtered_adj_phased_gntp

bcftools view \
    --include 'TYPE=="snp" && INFO/MAF>=0.0005 && INFO/R2>=0.3 && N_ALT==1' \
    --threads $runCpus \
    --output-type z \
    --output-file \$imputed_filtered_adj_phased_gntp \
    \${imputed_filtered_adj_phased_gntp/.vcf.gz/-singlealt.vcf.gz}

rm -f \$imputed_filtered_adj_phased_gntp
EOF

