#!/bin/bash
#
## WASP pipeline. Remapping bias reads
#

set -o errexit
set -o errtrace

# WASP snp2h5. Construct a HDF5 database from VCF file including haplotypes
# --geno_prob $geno_prob_file \   to use this option, the vcf file has to include GL or GP token

module purge
module load HDF5/1.8.14-foss-2015b
module list

mkdir -p $hdf5_dir
$snp2h5_exe \
	--chrom $chrm_info \
	--format vcf \
	--snp_tab $snp_tab_file \
	--haplotype $haplotype_file \
	--snp_index $snp_index_file \
	$smpl_snps
