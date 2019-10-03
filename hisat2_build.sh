#!/bin/bash
#SBATCH --time=2:59:50
#SBATCH --mem=180G
#SBATCH --cpus=15
#SBATCH --err=%j-%u-hisat_build.err
#SBATCH --output=%j-%u-hisat_build.log

module load Python/3.6.3-foss-2015b

# Project directories
my_proj_dir=~/Documents/projects/ASECausalSNPPrioritization
my_ipt_dir=$my_proj_dir/inputs
my_buf_dir=$my_proj_dir/buffers
my_tmp_dir=$my_proj_dir/temps

# The executable file
my_tlkt_dir=/groups/umcg-gcc/tmp03/umcg-zzhang/tools
my_hisat2_build_exe=$my_tlkt_dir/hisat2-2.1.0/hisat2-build
my_hess_py=$my_tlkt_dir/hisat2-2.1.0/hisat2_extract_splice_sites.py
my_hee_py=$my_tlkt_dir/hisat2-2.1.0/hisat2_extract_exons.py

# Input
my_gff_file=$my_ipt_dir/GRCh37_reference/Homo_sapiens.GRCh37.87.gtf
my_ref_genome=$my_ipt_dir/GRCh37_reference/human_g1k_v37.fasta

# Output
my_idx_dir=$my_buf_dir/hisat2_ss_exn_idx  # With splice sites and exons

# Extract splicing sites
$my_hess_py $my_gff_file > $my_tmp_dir/splices_sites.tsv
$my_hee_py $my_gff_file > $my_tmp_dir/exons.tsv

# Buid index
$my_hisat2_build_exe \
    -p 15 \
    --ss $my_tmp_dir/splices_sites.tsv \
    --exon $my_tmp_dir/exons.tsv \
    $my_ref_genome \
    ${my_ref_genome%.*}

rm $my_tmp_dir/splices_sites.tsv $my_tmp_dir/exons.tsv