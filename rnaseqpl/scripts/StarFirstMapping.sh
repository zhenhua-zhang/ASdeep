#!/bin/bash
#
## STAR alignReads. First round mapping FASTQ reads to reference genome.
#

# TODO:
#     1. A better allocation of threads for runThreadN and outBAMsortingThreadN
#     2. Specify genomeDir

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

Note:
    1. The name of input FASTQ files regex \${fastqPrefix}\${fastqId}_R{1,2}\${fastqSuffix}

More information please contact Zhenhua Zhang <zhenhua.zhang217@gmail.com>
EOF
    exit 0
}

opt=$(getopt -l "workDir:,fastqId:,fastqDir:,fastqPrefix:,fastqSuffix:,help" -- "w:i:g:p:P:S:" $@)
eval set -- ${opt}
while true; do
    case $1 in
        -w|--workDir) shift && workDir=$1 ;;
        -g|--genomeDir) shift && genomeDir=$1 ;;
        -i|--fastqId) shift && fastqId=$1 ;;
        -p|--fastqDir) shift && fastqDir=$1 ;;
        -S|--fastqSuffix) shift && fastqSuffix=$1 ;;
        -P|--fastqPrefix) shift && fastqPrefix=$1 ;;
        --outBAMsortingThreadN) shift && outBAMsortingThreadN=$1 ;;
        --help) echo_help ;;
        --) echo_help;;
    esac
    shift
done

module purge
module load STAR/2.5.1b-foss-2015b
module list

workDir=${workDir:?[E]: -w/--workDir is required!}
genomeDir=${genomeDir:?[E]: -g/--genomeDir is required!}

fastqId=${fastqId:?[E]: -i/--fastqId is required!}
fastqDir=${fastqDir:?[E]: -p/--fastqDir is required!}

readFilesIn_1=${fastqDir}/${fastqPrefix}${fastqId}"_R1"${fastqSuffix}
readFilesIn_2=${fastqDir}/${fastqPrefix}${fastqId}"_R2"${fastqSuffix}

outBAMsortingThreadN=${outBAMsortingThreadN:=2}
[ -n ${SLURM_CPUS_PER_TASK} ] \
    && runThreadN=$[ ${SLURM_CPUS_PER_TASK} - ${outBAMsortingThreadN} ]
    || runThreadN=$[ $(grep -c processor /proc/cpuinfo) - ${outBAMsortingThreadN} ]

outSAMattrRGline="ID:${fastqId} SM:${fastqId} PL:ILLUMINA"
outFileNamePrefix=${workDir}/tmpdir/${fastqId}/starTmpdir

STAR \
    --runMode alignReads \
    --genomeDir ${genomeDir} \
    --runThreadN ${runThreadN} \
    --readFilesIn ${readFilesIn_1} ${readFilesIn_2} \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattrRGline ${outSAMattrRGline} \
    --outBAMsortingThreadN ${outBAMSortingThreadN} \
    --twopassMode "Basic" \
    --outFileNamePrefix ${outFileNamePrefix}


module purge
module load SAMtools
module list

samtools index ${outFileNamePrefix}*.bam
