#!/bin/bash
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: 2021-04-24

#
## A utility script to preprocess variants data for GEUVADIS cohort
#

# Update the resource parameters

reqmem=5G
reqtime=02:59:00
reqcpus=5
output=
jobname=

sbatch --mem $reqmem \
  --time $reqtime \
  --cpus-per-task $reqcpus \
  --job-name geuvadis_genotype_qc \
  --output %j.geuvadis_genotype_qc.log \
  <<'EOF'
#!/bin/bash
source /apps/modules/modules.bashrc

# Should be after the loading of modules.bashrc.
set -Ee -o pipefail

proj_dir=$HOME/Documents/projects/wp_ase_dlp
smpfile=$proj_dir/inputs/Geuvadis/samples/igsr_Geuvadis-mRNA_Phase3_sample-name.txt
test ! -e $smpfile && echo "Failed to find file: line ${LINENO}" && exit -1

module purge
module load BCFtools
module list

in_dir=$proj_dir/inputs/Geuvadis/genotypes
out_dir=$proj_dir/outputs/haplotypes
for chr in {1..22}; do
    in_vcf=$in_dir/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
    out_vcf=$out_dir/Geuvadis.chr${chr}.vcf.gz

    bcftools view \
        --threads 5 \
        --output-type z \
        --min-alleles 2 \
        --max-alleles 2 \
        --force-samples \
        --output-file $out_vcf \
        --samples-file $smpfile \
        --types snps \
        --include 'INFO/EUR_AF >= 0.005' \
        $in_vcf

    bcftools index --threads 5 --tbi --force $out_vcf
done

# Using 1..22 to keep variants sorted by chromosome
bcftools concat \
    --allow-overlaps \
    --output $out_dir/all-Geuvadis-singleAlt.vcf.gz \
    --output-type z \
    --threads 5 \
    $out_dir/Geuvadis.chr{1..22}.vcf.gz

bcftools index -t $out_dir/all-Geuvadis-singleAlt.vcf.gz

rm -fr $out_dir/Geuvadis.*

echo Job done!
EOF
