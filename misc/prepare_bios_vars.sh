#!/bin/bash
#
# File name : prepare_bios_vars.sh
# Author    : Zhenhua Zhang
# E-mail    : zhenhua.zhang217@gmail.com
# Created   : Tue 30 Jun 2020 04:49:32 PM CEST
# Version   : V 0.1.0
# License   : MIT
#

# TODO: Recover to the first version.

if [[ -e /apps/modules/modules.bashrc ]]; then
    source /apps/modules/modules.bashrc
fi

set -Eeu -o pipefail

debug=true
if [[ $debug == true ]]; then
    runCpus=2
    runMems=4G
    runTime=0:59:0
    runArray=22-22
else
    runCpus=2
    runMems=15G
    runTime=19:59:0
    runArray=1-22
fi

sbatch \
    --mem=$runMems \
    --time=$runTime \
    --array=$runArray \
    --cpus-per-task=$runCpus \
    --job-name=phasing \
    --output=%A-%a_%u_phasing.log \
    <<"EOF" | cut -f4 -d' '
#!/bin/bash

# Ensure Slurm settings is updated.
if [[ -e /apps/modules/modules.bashrc ]]; then
    source /apps/modules/modules.bashrc
fi

# Set some shell attributes
set -Eeu -o pipefail

# Project (working) dirs
pjdir=/groups/umcg-bios/tmp01/users/umcg-zzhang/projects/wp_ase_dlp

# Version of phasing
phasing_version=haplotypes

# Working on per chrommsome by jobarray of SLURM
chrid=${SLURM_ARRAY_TASK_ID:-22}

# Threads environmental variable from SLURM
# Multiple CPUs per task, but --cpus-per-task is required to use the SLURM_CPUS_PER_TASK
# variable
threads=${SLURM_CPUS_PER_TASK:-$(grep -c processor /proc/cpuinfo)}

# Files for each step
bcftools_dir=$pjdir/outputs/$phasing_version/$chrid/bcftools
conformGt_dir=$pjdir/outputs/$phasing_version/$chrid/conformGt
beagle_dir=$pjdir/outputs/$phasing_version/$chrid/beagle

# Create buffer directories
mkdir -p $bcftools_dir $conformGt_dir $beagle_dir

# BCFtools to filter the variants
module purge && module load BCFtools/1.9-foss-2018b && module list

# Input: upstream
# The imputed files are located (calculon): /groups/umcg-bios/prm02/projects/HRC_imputation,
# which are imputed using Michigan Imputation Server.
id_mapping_file=$pjdir/inputs/idmapping/BIOS-genotype-rnaseq-ids_with-LL_usable_20210105.txt

# Output: downstream
imputed_filtered_gntp=$bcftools_dir/bios-${chrid}.filtered.merged.vcf.gz

for vcffile in $pjdir/inputs/BIOS_genotypes/*/chr${chrid}.dose.vcf.gz; do
    cohort=$(basename $(dirname $vcffile))
    vcffile_bnm=$cohort-$(basename $vcffile)

    # Rename samples: using RNA-seq id instead of genotyping id
    vcffile_rh=$bcftools_dir/${vcffile_bnm/.vcf.gz/.rh.vcf.gz}
    bcftools reheader \
        --samples <(cut -f2,3 $id_mapping_file) \
        --output $vcffile_rh \
        $vcffile

    bcftools index \
        --threads $threads \
        --force \
        --tbi \
        $vcffile_rh

    # Select only samples in id mapping file, i.e., samples with both RNA-seq and genotypes
    # R2 refer to Xiaojing's QC on 200HIV genotypes (>=0.3);
    vcffile_rh_ni=$bcftools_dir/${vcffile_bnm/.vcf.gz/.rh.ni.vcf.gz}
    bcftools view \
        --force-samples \
        --samples-file <(cut -f3 $id_mapping_file) \
        --output-type z \
        --output-file $vcffile_rh_ni \
        --threads $threads \
        $vcffile_rh
        # --include 'INFO/R2 >= 0.3' \

    # Index it as bcftools/merge requires index
    bcftools index \
        --threads $threads \
        --force \
        --tbi \
        $vcffile_rh_ni

    rm -f $vcffile_rh $vcffile_rh.tbi
done

# Merge all usable variants. 
# NOTE: the INFO/MAF should be updated and INFO[R2,[ER2]] fields should be ignored
# TODO: A way to allocate CPUs for the used subcmds including merge, annotate, plugin, view.
bcftools merge \
    --force-samples \
    --threads $(( $threads - 3 <= 0 ? 0 : $threads - 3 )) \
    $bcftools_dir/*.rh.ni.vcf.gz \
    | bcftools annotate --remove INFO/MAF,INFO/R2,INFO/ER2 \
    | bcftools plugin fill-tags -- -t MAF,AF \
    | bcftools view \
        --include 'TYPE=="snp" && INFO/MAF >= 0.005 && N_ALT == 1' \
        --output-type z \
        --output-file $imputed_filtered_gntp \

# Beagle version 5.1 which is powered by Java 1.8 or higher.
module purge && module load Java/1.8.0_144 && module list

# Input: reference
ref_panel_gntp=$pjdir/inputs/beagle/reference_panel/chr${chrid}.1kg.phase3.v5a.vcf.gz

# Input: upstream
# imputed_filtered_gntp

# Output: downstream
imputed_filtered_adj_gntp=$conformGt_dir/$(basename ${imputed_filtered_gntp/vcf.gz/adjByRef})

java -jar ~/tools/beagle/conform-gt.24May16.cee.jar \
    gt=$imputed_filtered_gntp \
    ref=$ref_panel_gntp \
    chrom=$chrid \
    match=POS \
    out=$imputed_filtered_adj_gntp

# Phasing.
# Input: reference map and genotype.
ref_panel_map=$pjdir/inputs/beagle/genetic_maps/plink.chr${chrid}.GRCh37.map

# Output: final
imputed_filtered_adj_phased_gntp=$beagle_dir/beagle-chr${chrid}

java -jar ~/tools/beagle/beagle.18May20.d20.jar \
    gt=$imputed_filtered_adj_gntp".vcf.gz" \
    ref=$ref_panel_gntp \
    map=$ref_panel_map \
    impute=false \
    seed=31415 \
    nthreads=$threads \
    out=$imputed_filtered_adj_phased_gntp

rm -f $imputed_filtered_adj_phased_gntp
EOF

