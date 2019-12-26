#!/bin/bash
#
## Preprocessing reads using fastp
#

# TODO:
#     1. Command lind processing script scope

# NOTE:
# --detect_adapter_for_pe is enabled to check the adapters, it could slow down the process by the manual.
# --trim_poly_g is enabled to trim poly-G, where the G means no signal in the Illumina two-color systems.
# --overrepresentation_sampling is used to check 1% of all reads to analysis the overrepresented sequence

source ../bin/utils

echo_help() {
    cat <<EOF

Help:
  -w|--workDir Required
      Working directory
  -i|--fastqId Required
      FASTQ files ID
  -p|--fastqDir Required
      Path to the FASTQ files
  -s|--fastqSuffix Optional
      Suffix for the input FASTQ files
  --threads Optional
      Number of threadings for fastp
  --help Optional
      Print this help context and exit.

More information please contact Zhenhua Zhang <zhenhua.zhang217@gmail.com>
EOF
    exit 0
}

opt=$(getopt -l "workDir:,fastqId:,fastqDir:,fastqSuffix:,help" -- "w:i:p:s:" $@)
eval set -- ${opt}
while true; do
    case $1 in
        -w|--workDir) shift && workdir=$1 ;;
        -i|--fastqId) shift && fastqId=$1 ;;
        -p|--fastqDir) shift && fastqDir=$1 ;;
        -s|--fastqSuffix) shift && fastqSuffix=$1 ;;
        --threads) shift && fastpThreads=$1 ;;
        --help) echo_help ;;
        --) echo_help;;
    esac
    shift
done

# Variables by command line
workDir=${workDir:?[E]: -w/--workDir is mandatory}
fastqId=${fastqId:?[E]: -i/--fastqId is mandatory}
fastqDir=${fastqDir:?[E]: -p/--fastqDir is mandatory}
fastqPrefix=${fastqPrefix:=}
fastqSuffix=${fastqSuffix:=.fq.gz}

# Variables by hard coding
[ -n $SLURM_CPUS_PER_TASK ] \
    && fastpThreads=$[ $SLURM_CPUS_PER_TASK ]
    || fastpThreads=$[ $(cat /proc/cpuinfo | grep -c processor) ]
fastpTmpDir=$workDir/$fastqId/fastpTmpDir

# Inputs
fastqFile_1=${fastqDir}/${fastqPrefix}${fastqId}"_R1"${fastqSuffix}
fastqFile_2=${fastqDir}/${fastqPrefix}${fastqId}"_R2"${fastqSuffix}

# Outputs
fastqPaired_1=$fastpTmpDir/${fastqId}_paired_R1.fq.gz
fastqPaired_2=$fastpTmpDir/${fastqId}_paired_R2.fq.gz

fastqHtmlReport=$fastpTmpDir/${fastqId}_report.html
fastqJsonReport=$fastpTmpDir/${fastqId}_report.json

fastqFailed=$fastpTmpDir/${fastqId}_failed.fq.gz
fastqUnpaired_1=$fastpTmpDir/${fastqId}_unpaired_R1.fq.gz
fastqUnpaired_2=$fastpTmpDir/${fastqId}_unpaired_R2.fq.gz

# Command line
../bin/fastp/fastp \
    --in1 $fastqFile_1 \
    --in2 $fastqFile_2 \
    --out1 $fastqPaired_1 \
    --out2 $fastqPaired_2 \
    --unpaired1 $fastqUnpaired_1 \
    --unpaired2 $fastqUnpaired_2  \
    --failed_out $fastqFailed \
    --html $fastqHtmlReport \
    --json $fastqJsonReport \
    --detect_adapter_for_pe \
    --cut_front \
    --cut_tail \
    --correction \
    --trim_poly_g \
    --trim_poly_x \
    --overrepresentation_sampling 100 \
    --thread $fastpThreads \
    --dont_overwrite
