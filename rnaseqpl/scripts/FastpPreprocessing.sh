#!/bin/bash
#
## Preprocessing reads using fastp
#

# NOTE:
# --detect_adapter_for_pe is enabled to check the adapters, it could slow down the process by the manual.
# --trim_poly_g is enabled to trim poly-G, where the G means no signal in the Illumina two-color systems.
# --overrepresentation_sampling is used to check 1% of all reads to analysis the overrepresented sequence

set -o errexit
set -o errtrace

source /apps/modules/modules.bashrc
source ${WASPPL_SCRIPT_PATH}/utils.sh

echo_help() {
    cat <<EOF

Help:
  -w|--workDir [Required]
      Working directory.
  -i|--fastqId [Required]
      FASTQ files ID.
  -p|--fastqDir [Required]
      Path to the FASTQ files.
  -s|--fastqSuffix [Optional]
      Suffix for the input FASTQ files. Default: .fq.gz
  -P|--fastpExe [Optional]
      The path to the executbale fastp. Default: ~/tools/fastp/fastp
  --threads [Optional]
      Number of threadings for fastp. Default: all
  --help
      Print this help context and exit.

More information please contact Zhenhua Zhang <zhenhua.zhang217@gmail.com>

EOF
    exit 0
}

long_opts="workDir:,fastqId:,fastqDir:,fastqSuffix:,fastpExe:,threads:,help" 
opt=$(getopt -l $long_opts -- "w:i:q:s:P:" $@)
eval set -- ${opt}
while true; do
    case $1 in
        -w|--workDir) shift && workDir=$1 ;;
        -i|--fastqId) shift && fastqId=$1 ;;
        -P|--fastpExe) shift && fastpExe=$1 ;;
        -q|--fastqDir) shift && fastqDir=$1 ;;
        -s|--fastqSuffix) shift && fastqSuffix=$1 ;;
        --threads) shift && fastpThreads=$1 ;;
        --help) echo_help ;;
        --) shift && break ;;
    esac
    shift
done

# Variables by command line
workDir=${workDir:?[E]: -w/--workDir is mandatory}
fastqId=${fastqId:?[E]: -i/--fastqId is mandatory}
fastqDir=${fastqDir:?[E]: -p/--fastqDir is mandatory}
fastpExe=${fastpExe=~/tools/bin/fastp}
fastqPrefix=${fastqPrefix:=}
fastqSuffix=${fastqSuffix:=.fq.gz}

fastpTmpDir=$workDir/tmpdir/$fastqId/fastpTmpdir
fastpThreads=$(check_cpus)

# Inputs
fastqFile_1=${fastqDir}/${fastqPrefix}${fastqId}"_R1"${fastqSuffix}
fastqFile_2=${fastqDir}/${fastqPrefix}${fastqId}"_R2"${fastqSuffix}

# Outputs
fastqPaired_1=$fastpTmpDir/${fastqId}_R1_paired.fq.gz
fastqPaired_2=$fastpTmpDir/${fastqId}_R2_paired.fq.gz

fastqHtmlReport=$fastpTmpDir/${fastqId}_report.html
fastqJsonReport=$fastpTmpDir/${fastqId}_report.json

fastqFailed=$fastpTmpDir/${fastqId}_failed.fq.gz
fastqUnpaired_1=$fastpTmpDir/${fastqId}_R1_unpaired.fq.gz
fastqUnpaired_2=$fastpTmpDir/${fastqId}_R2_unpaired.fq.gz

# Command line
$fastpExe \
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

