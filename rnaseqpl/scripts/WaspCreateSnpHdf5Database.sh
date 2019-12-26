#!/bin/bash
#
## WASP pipeline. Remapping bias reads
#

set -o errexit
set -o errtrace

source ../bin/utils

echo_version() {
    cat << EOF

$(basename $0), Version ${SCRIPT_VERSOIN:=UNKNOWN}
EOF
}

echo_usage() {
    cat <<EOF

Usage: ./$0 [options]
    or bash $0 [options]
EOF
}

echo_help() {
    cat <<EOF
$(basename $0), Version ${SCRIPT_VERSOIN:=UNKNOWN}

Help:
  -h, --help    Optional. Action: print_info
    Print this help context and exit.
  -u, --usage    Optional. Action: print_info
    Print usage context and exit.
  -V, --version    Optional. Action: store_true
    Print version of current script and exit.

More information please contact Zhenhua Zhang <zhenhua.zhang217@gmail.com>
EOF
}

opt=$(getopt -l "workDir:,fastqId:,snpHdf5db:,help,usage,version" -- "w:i:p:g:f:a:v:" $@)
eval set -- ${opt}
while true; do
    case $1 in
        -w|--workDir) shift && workdir=$1 ;;
        -v|--vcfFile) shift && vcfFile=$1 ;;  # the SNP hDF5 database
        -c|--chromInfoFile) shift && chrmInfoFile=$1 ;;  # the SNP hDF5 database
        -e|--snp2h5Exe) shift && snp2h5Exe=$1 ;;
        --help) echo_help && exit 0;;
        --usage) echo_usage && exit 0;;
        --version) echo_version && exit 0;;
        --) echo_help && exit 0;;
        *) ehco_help && exit 0;;
    esac
    shift
done

# WASP snp2h5. Construct a HDF5 database from VCF file including haplotypes
# --geno_prob $geno_prob_file \   to use this option, the vcf file has to include GL or GP token

workDir=${workDir:?-w/--workDir is required}
fastqId=${fastqId:?-s/--fastqId is required}  # ${fastqId}_R1.fq.gz & ${fastqId}_R2.fq.gz

module purge
module load HDF5/1.8.14-foss-2015b
module list

# Create a directory for SNP2H5 and `cd` into it.
snpHDF5db=${workDir}/snpHDF5db
mkdir -p ${snpHDF5db}

${snp2h5Exe:=snp2h5} \
	--chrom ${chrmInfoFile} \
	--format vcf \
	--snp_tab ${snpHDF5db}/snpsTab.h5 \
	--snp_index ${snpHDF5db}/snpsIndex.h5 \
	--haplotype ${snpHDF5db}/haplotype.h5 \
	${vcfFile}

[ $? -eq 0 ] && INFO "Job was done!" || ERRO "Job exit with non-zero!"
