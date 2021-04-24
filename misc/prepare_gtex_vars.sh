#!/bin/bash
# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: 2021-04-24

#
## A utility script to preprocessing variats data for GTEx cohort
#

# Project directories
projdir=~/Documents/projects/wp_ase_dlp
iptdir=$projdir/inputs
optdir=$projdir/outputs

# Filter variants
raw_vars=$iptdir/GTEx/vcf/phASER_WASP_GTEx_v8_merged.vcf.gz
clean_vars=${raw_vars/.vcf.gz/.clean.nochr.vcf}

# Liftover
clean_vars_b37=${clean_vars/.vcf.gz/.b37.vcf}
ref_genome_b37=$iptdir/GRCh37_reference/human_g1k_v37.fasta
convert_chain=$iptdir/Ensembl_references/GRCh38_to_GRCh37.chain.gz

# Update the resource parameters.
stepmem=1G
stepcpus=5
steptime=2:59:59

sbatch \
    --mem $stepmem \
    --time $steptime \
    --cpus-per-task $stepcpus \
    --job-name filter-gtex-vars \
    --output %j-filter_gtex_vars.log \
    <<EOF
#!/bin/bash
set -Ee -o pipefail
module load BCFtools/1.9-foss-2018a

# Generate index
bcftools index \
    --threads $stepcpus $raw_vars

# Filter variants
bcftools view 
    --threads $stepcpus \
    -m 2 \
    -M 2 \
    -v snps \
    -i 'INFO/AF > 0.005 && INFO/AF <= 0.995' \
    $raw_vars $(echo chr{1..22} | tr ' ' ',') \
    | bcftools annotate --rename-chrs $chr_mapping_file \
    >| ${clean_vars}

# Liftover from GRCh38 to GRCh37
module purge
module load Python/3.7.4-GCCcore-8.3.0
source $projdir/script/.env/bin/activate
CrossMap.py vcf $convert_chain $clean_vars $ref_genome_b37 $clean_vars_b37 --no-comp-alleles

# Compress VCF
module purge
module load bcftools/1.9-foss-2018a
bgzip -f -i -@ $stepcpus $clean_vars_b37
bcftools index -@ $stepcpus -t $clean_vars_b37

# Clean up
rm -f ${clean_vars}

# Report
echo Done
EOF

