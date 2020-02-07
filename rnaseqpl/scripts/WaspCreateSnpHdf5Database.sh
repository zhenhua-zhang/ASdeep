#!/bin/bash
#
## WASP pipeline. Remapping bias reads
#

# TODO:
#   1. Add --samples options to specify samples

set -o errexit
set -o errtrace

source /apps/modules/modules.bashrc

ERRO() {
    echo -e "[E]: $1" >&2 && exit -1
}

WARN() {
    echo -e "[W]: $1" >&2
}

INFO() {
    echo -e "[I]: $1"
}

check_cpus() {
    local n_cpus
    local balance

    [ $1"x" == "x" ] \
        && balance=0 \
        || balance=$1

    [ $SLURM_CPUS_PER_TASK"x" == "x" ] \
        && n_cpus=$[ $(grep -c processor /proc/cpuinfo) - $balance ] \
        || n_cpus=$[ $SLURM_CPUS_PER_TASK - $balance ]

    [ $n_cpus -gt 0 ] \
        && echo $n_cpus \
        || echo $[ $n_cpus + $balance ]
}

echo_help() {
    cat <<EOF

Usage: $(basename $0) [options]

Help:
  -w|--workDir [required]
      Working directory.
  -c|--chromInfoFile [required]
      Length of each chromosome.
  -e|--snp2h5Exe [Optional]
      The executbale snp2h5 programe. For currently the WASP pipeline hasn't
      globally deployed.
  -h, --help
    Print this help context and exit.

More information please contact Zhenhua Zhang <zhenhua.zhang217@gmail.com>

EOF
}

opt=$(getopt -l "workDir:,vcfFileDir:,vcfFilePre:,vcfFileSuf:,chromInfoFile:,snp2h5Exe:,help" -- "w:v:c:e:h" $@)
eval set -- ${opt}
while true; do
    case $1 in
        -w|--workDir) shift && workDir=$1 ;;
        -v|--vcfFileDir) shift && vcfFileDir=$1 ;;
        -P|--vcfFilePre) shift && vcfFilePre=$1 ;;
        -S|--vcfFileSuf) shift && vcfFileSuf=$1 ;;
        -c|--chromInfoFile) shift && chromInfoFile=$1 ;;
        -e|--snp2h5Exe) shift && snp2h5Exe=$1 ;;
        -s|--samples) shift && samples=$1 ;;
        -h|--help) echo_help && exit 0;;
        --) shift && break ;;
    esac
    shift
done

# WASP snp2h5. Construct a HDF5 database from VCF file including haplotypes
# --geno_prob $geno_prob_file \   to use this option, the vcf file has to include GL or GP token
chromId=${SLURM_ARRAY_TASK_ID:=22}
workDir=${workDir:?-w/--workDir is required!}
vcfFileDir=${vcfFileDir:?-v/--vcfFileDir is required!}
chromInfoFile=${chromInfoFile:?-c/--chromInfoFile is required!}
vcfFilePre=${vcfFilePre:=gonl.chr}
vcfFileSuf=${vcfFileSuf:=.snps.r5.3.vcf.gz}

vcfFile=$vcfFileDir/${vcfFilePre}${chromId}${vcfFileSuf}

module purge
module load HDF5/1.8.14-foss-2015b
module list

# Create a directory for SNP2H5
snph5db=$workDir/snph5db
mkdir -p $snph5db/$chromId

${snp2h5Exe:=snp2h5} \
	--chrom $chromInfoFile \
	--format vcf \
	--snp_tab $snph5db/$chromId/snps_tab.h5 \
	--snp_index $snph5db/$chromId/snps_index.h5 \
	--haplotype $snph5db/$chromId/haplotype.h5 \
	$vcfFile
