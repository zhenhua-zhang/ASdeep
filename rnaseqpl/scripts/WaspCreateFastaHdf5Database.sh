#!/bin/bash
#
## WASP pipeline. Remapping bias reads
#

set -o errexit
set -o errtrace

source /apps/modules/modules.bashrc
source ${WASPPL_SCRIPT_PATH}/utils.sh

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

long_opts="workDir:,sequenceFile:,chromInfoFile:,fasta2h5Exe:,help"
opt=$(getopt -l $long_opts -- "w:s:c:e:h" $@)
eval set -- ${opt}
while true; do
    case $1 in
        -w|--workDir) shift && workDir=$1 ;;
        -s|--sequenceFile) shift && sequenceFile=$1 ;;
        -c|--chromInfoFile) shift && chromInfoFile=$1 ;;
        -e|--fasta2h5Exe) shift && fasta2h5Exe=$1 ;;
        -h|--help) echo_help && exit 0;;
        --) shift && break ;;
    esac
    shift
done

# WASP fasta2h5. Construct a HDF5 database from genome sequence
workDir=${workDir:?-w/--workDir is required!}
sequenceFile=${sequenceFile:?-s/--sequenceFile is required!}
chromInfoFile=${chromInfoFile:?-c/--chromInfoFile is required!}

module purge
module load HDF5/1.8.14-foss-2015b
module list

# Create a directory for fasta
fastah5db=$workDir/fastah5db
mkdir -p $fastah5db

csplit -s -f $fastah5db/chr -b '%d' ${sequenceFile} '/>/' {22} && rm -f chr{0,23}

for chr_file in $(ls $fastah5db/chr*); do
    ${fasta2h5Exe:=fasta2h5} \
        --chrom $chromInfoFile \
        --seq $fastah5db/${chr_file}.h5 \
        $fastah5db/${chr_file}
    rm ${chr_file} -f
done
