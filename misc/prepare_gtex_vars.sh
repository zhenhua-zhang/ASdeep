#!/bin/bash
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: 2021-04-24

#
## A utility script to preprocessing variats data for GTEx cohort
#

# Update the resource parameters.
reqmem=32G
reqcpus=1
reqtime=7:59:59

sbatch \
    --mem $reqmem \
    --time $reqtime \
    --cpus-per-task $reqcpus \
    --job-name filter-gtex-vars \
    --output %j-filter_gtex_vars.log \
    <<'EOF'
#!/bin/bash
set -Ee -o pipefail

#
## Project directories, computational resource, input files
#
# Project directories
projdir=~/Documents/projects/wp_ase_dlp
iptdir=$projdir/inputs
optdir=$projdir/outputs

# Computational resource
threads=${SLURM_CPUS_PER_TASK:-$(grep -c processor /proc/cpuinfo)}

# Raw variants from GTEx
raw_vars=$iptdir/GTEx/vcf/phASER_WASP_GTEx_v8_merged.vcf.gz


#
## Generate index
#
module purge && module load BCFtools/1.9-foss-2018a && module list

if [[ ! -f $raw_vars.tbi ]]; then
    bcftools index --threads $threads --tbi $threads $raw_vars
else
    echo Find index for $raw_vars
fi


#
## Filter variants. hq = high quality; nc = no chr
#
# Filter variants. Only biallelic SNPs with AF from 0.005 to 0.995 and only 1 alternative allele
chosen_hq_vars=$optdir/phasing/gtex/$(sed 's/.vcf.gz/_hq.vcf/g' <<<$(basename $raw_vars))
chr_mapping_file=$iptdir/GTEx/vcf/phASER_WASP_GTEx_v8_merged.chr.txt

bcftools view \
    -m 2 \
    -M 2 \
    -i 'TYPE == "snp" && INFO/AF > 0.005 && INFO/AF <= 0.995 && N_ALT == 1' \
    --threads $(( $threads / 2 )) \
    $raw_vars $(echo chr{1..22} | tr " " ",") \
    | bcftools annotate \
        -O v \
        -o $chosen_hq_vars \
        --threads $(( $threads / 2 )) \
        --rename-chrs $chr_mapping_file


#
## Liftover from GRCh38 to GRCh37
#
module purge && module load picard && module list

# Inputs and outputs for picard/LiftoverVcf
chosen_hq_nc_b37_vars=${chosen_hq_vars/.vcf/_b37.vcf}
reject_vars=${chosen_hq_vars/.vcf/_rej.vcf}
ref_genome_b37=$iptdir/GRCh37_reference/human_g1k_v37.fasta
convert_chain=$iptdir/UCSC_browser/hg38ToHg19.over.nochr.chain.gz

java -jar $EBROOTPICARD/picard.jar LiftoverVcf \
    -I $chosen_hq_vars \
    -O $chosen_hq_nc_b37_vars \
    -R $ref_genome_b37 \
    -CHAIN $convert_chain \
    -REJECT $reject_vars \
    -WARN_ON_MISSING_CONTIG true \
    --TMP_DIR $optdir/phasing/gtex


#
## Annotate VCF
#
module purge && module load BCFtools/1.9-foss-2018a && module list

# Input and output
ref_vars=$iptdir/dbSNP/00-All.vcf.gz
chosen_hq_nc_b37_annot_vars=${chosen_hq_nc_b37_vars}_annot.vcf.gz

# Compress successfully liftovered varaints and generate index
bcftools annotate \
    -O z \
    -c ID \
    -a $ref_vars \
    --threads $threads \
    -o $chosen_hq_nc_b37_annot_vars \
    $chosen_hq_nc_b37_vars

bcftools index --threads $threads --tbi $chosen_hq_nc_b37_annot_vars

# Compress rejected variants
bcftools view \
    --threads $threads \
    --output-type z \
    --output-file $reject_vars.gz \
    $reject_vars


#
## Clean up
#
rm -f $chosen_hq_vars $chosen_hq_nc_b37_vars $chosen_hq_nc_b37_vars.idx $reject_vars
echo Done

EOF

