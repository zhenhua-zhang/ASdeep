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

if [[ -e /apps/modules/modules.bashrc ]]; then
    source /apps/modules/modules.bashrc
fi

# Cohorts: CODAM  NTR  PAN  RS
# Sample information
cohort="${1:?Error: You should give the cohort id!}"

# Sub dirs
#pjdir=/Documents/projects/wp_ase_dlp
pjdir=/groups/umcg-bios/tmp03/users/umcg-zzhang/projects/wp_ase_dlp
if [[ ! -d $pjdir/inputs/BIOS_genotypes/$cohort ]]; then
    echo "Failed to find $pjdir/inputs/BIOS_genotypes/$cohort! Exit"
    exit
fi

bcftools_dir=$pjdir/outputs/phasing/$cohort/bcftools
conformGt_dir=$pjdir/outputs/phasing/$cohort/conformGt
beagle_dir=$pjdir/outputs/phasing/$cohort/beagle
log_dir=$pjdir/outputs/phasing/$cohort/logs
mkdir -p $bcftools_dir $conformGt_dir $beagle_dir $log_dir

debug=false
if [[ $debug == true ]]; then
    runCpus=2
    runMems=15G
    runTime=0:29:0
    runArray=22-22
    set_flags=-Eeux
else
    runCpus=15
    runMems=40G
    runTime=19:59:0
    runArray=1-22
    set_flags=-Eeu
fi

# Beagle settings
seed=31415


sbatch --time=$runTime \
    --mem=$runMems \
    --output=$log_dir/%A-%a_%u_phasing_$cohort.log \
    --job-name=phasing-$cohort \
    --array=$runArray \
    --cpus-per-task=$runCpus \
    <<EOF | cut -f4 -d' '
#!/bin/bash
# Ensure Slurm settings is updated.
source /apps/modules/modules.bashrc

set $set_flags -o pipefail

chr_id=\${SLURM_ARRAY_TASK_ID:=1}

cat <<INNEREOF
#
## Quality control by BCFtools
#
INNEREOF
module purge && module load BCFtools && module list

# BCFtools version should be: bcftools=1.7, htslib=1.7
bcftools --version

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


cat <<INNEREOF
#
## Phasing by Beagle
#
INNEREOF

# Beagle version 5.1 which is powered by Java 1.8 or higher.
module purge && module load Java/1.8.0_45 && module list

# conform version;
java -jar ~/tools/beagle/conform-gt.24May16.cee.jar | grep usage

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

# Beagle version: beagle.18May20.d20.jar (version 5.1)
java -jar ~/tools/beagle/beagle.18May20.d20.jar | head -1

# Imputation and phasing.
# Input: reference map and genotype.
# ref_panel_gntp=$pjdir/inputs/beagle/reference_panel/chr\${chr_id}.1kg.phase3.v5a.vcf.gz
ref_panel_map=$pjdir/inputs/beagle/genetic_maps/plink.chr\${chr_id}.GRCh37.map

# Input: upstream
# imputed_filtered_adj_gntp".vcf.gz"

# Output: final
imputed_filtered_adj_phased_gntp=$beagle_dir/beagle-$cohort-chr\${chr_id}

java -jar ~/tools/beagle/beagle.18May20.d20.jar \
    gt=\$imputed_filtered_adj_gntp".vcf.gz" \
    ref=\$ref_panel_gntp \
    map=\$ref_panel_map \
    impute=false \
    seed=$seed \
    nthreads=$runCpus \
    out=\$imputed_filtered_adj_phased_gntp

EOF


# FIXME:
# 
# runTime=1:59:0
# runMems=40G
# log=$log_dir/%A_%a-%u-genotypePooling-$cohort.log
# job=phasingAndImputation-$cohort
# sbatch --time=$runTime \
#     --mem=$runMems \
#     --output=$log \
#     --job-name=$job \
#     --array=1-7 \
#     --cpus-per-task=${runCpus} \
#     <<EOF | cut -f4 -d' '
# #!/bin/bash
# if [[ -e /apps/modules/modules.bashrc ]]; then
#     source /apps/modules/modules.bashrc
# fi
# 
# set -Eeux -o pipefail
# 
# echo -e '#\n## Rename samples and merge VCF files\n#'
# 
# module purge && module load BCFtools && module list
# bcftools --version | head -2
# 
# echo -e '#\n## PCA analysis to identify population stratification\n#'
# module purge && module load plink && module list
# plink --version
# 
# EOF
