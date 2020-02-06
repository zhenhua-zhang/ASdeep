#!/bin/bash
#
## WASP pipeline. Remapping bias reads
#


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
  -e|--fasta2h5Exe [Optional]
      The executbale fasta2h5 programe. For currently the WASP pipeline hasn't
      globally deployed.
  -h, --help
    Print this help context and exit.

More information please contact Zhenhua Zhang <zhenhua.zhang217@gmail.com>

EOF
}

opt=$(getopt -l "workDir:,sequenceDir:,chromInfoFile:,fasta2h5Exe:,help" -- "w:s:c:e:h" $@)
eval set -- ${opt}
while true; do
    case $1 in
        -w|--workDir) shift && workDir=$1 ;;
        -s|--sequenceDir) shift && sequenceDir=$1 ;;
        -c|--chromInfoFile) shift && chromInfoFile=$1 ;;
        -e|--fasta2h5Exe) shift && fasta2h5Exe=$1 ;;
        -h|--help) echo_help && exit 0;;
        --) shift && break ;;
    esac
    shift
done

# WASP fasta2h5. Construct a HDF5 database from genome sequence
workDir=${workDir:?-w/--workDir is required!}
sequenceDir=${sequenceDir:?-s/--sequenceDir is required!}
chromInfoFile=${chromInfoFile:?-c/--chromInfoFile is required!}
chromId=${SLURM_ARRAY_TASK_ID:=22}

module purge
module load HDF5/1.8.14-foss-2015b
module list

# Create a directory for fasta
fastah5db=$workDir/fastah5db
mkdir -p $fastah5db

${fasta2h5Exe:=fasta2h5} \
	--chrom $chromInfoFile \
	--seq $fastah5db/chr${chromId}.h5 \
	$sequenceDir/chr${chromId}
