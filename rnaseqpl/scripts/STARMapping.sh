#!/bin/bash
#
## STAR alignReads. First round mapping FASTQ reads to reference genome.
#

# TODO:
#     1. A better allocation of threads for runThreadN and outBAMsortingThreadN
#     2. Specify genomeDir

# Note:
#     1. For human genome, it requires at least 35G memory.

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

Usage:
  bash $(basename $0) -w work_dir -g genome_dir -i fastq_id -p fastq_dir \\
    [-P fastq_prefix] [-S fastq_suffix] [--outBAMsortingThreadN n]

Help:
  -w|--workDir [Required]
      Working directory.
  -g|--genomeDir) [Required]
      Genome index path.
  -i|--fastqId [Required]
      FASTQ files ID.
  -p|--fastqDir [Required]
      Path to the FASTQ files.
  -o|--outFileNamePrefix [Optional]
      The prefix of the output file
  -P|--fastqPrefix [Optional]
      Prefix for the input FASTQ files.
  -s|--fastqSuffix [Optional]
      Suffix for the input FASTQ files.
  --outBAMsortingThreadN [Optional]
      Number of threadings to sort the BAM file.
  --help
      Print this help context and exit.

Note:
    1. The name of input FASTQ files regex \${fastqPrefix}\${fastqId}_R{1,2}\${fastqSuffix}

More information please contact Zhenhua Zhang <zhenhua.zhang217@gmail.com>

EOF
    exit 0
}

opt=$(getopt -l "workDir:,genomeDir:,fastqId:,fastqDir:,outFileNamePrefix:,fastqPrefix:,fastqSuffix:,fq1Pattern:,fq2Pattern:,outBAMsortingThreadN:,help" -- "w:g:i:p:o:P:S:" $@)
eval set -- ${opt}

while true; do
    case $1 in
        -w|--workDir) shift && workDir=$1 ;;
        -g|--genomeDir) shift && genomeDir=$1 ;;
        -i|--fastqId) shift && fastqId=$1 ;;
        -p|--fastqDir) shift && fastqDir=$1 ;;
        -o|--outFileNamePrefix) shift && outFileNamePrefix=$1 ;;
        -P|--fastqPrefix) shift && fastqPrefix=$1 ;;
        -S|--fastqSuffix) shift && fastqSuffix=$1 ;;
        --fq1Pattern) shift && fq1Pattern=$1 ;;
        --fq2Pattern) shift && fq2Pattern=$1 ;;
        --outBAMsortingThreadN) shift && outBAMsortingThreadN=$1 ;;
        --help) echo_help ;;
        --) shift && break ;;
    esac
    shift
done

workDir=${workDir:?[E]: -w/--workDir is required!}
workDir=$(readlink -f $workDir)

genomeDir=${genomeDir:?[E]: -g/--genomeDir is required!}
genomeDir=$(readlink -f $genomeDir)

fastqId=${fastqId:?[E]: -i/--fastqId is required!}

fastqDir=${fastqDir:?[E]: -p/--fastqDir is required!}
fastqDir=$(readlink -f $fastqDir)

fastqPrefix=${fastqPrefix:=}
fastqSuffix=${fastqSuffix:=.fq.gz}

if [ $fq1Pattern"xxx" == "xxx" ]; then
    readFilesIn_1=$fastqDir/$fastqPrefix$fastqId"_R1"$fastqSuffix
else
    readFilesIn_1=$(ls -m $fastqDir/$fq1Pattern | tr -d "\n ")
fi

if [ $fq2Pattern"xxx" == "xxx" ]; then
    readFilesIn_2=$fastqDir/$fastqPrefix$fastqId"_R2"$fastqSuffix
else
    readFilesIn_2=$(ls -m $fastqDir/$fq2Pattern | tr -d "\n ")
fi

outBAMsortingThreadN=${outBAMsortingThreadN:=2}
runThreadN=$(check_cpus $outBAMsortingThreadN)

outSAMattrRGline="ID:$fastqId SM:$fastqId PL:ILLUMINA"
outFileNamePrefix=${outFileNamePrefix:=$workDir/tmpdir/$fastqId/starTmpdir/$fastqId}

module purge
module load STAR/2.5.1b-foss-2015b
module list

STAR \
    --runMode alignReads \
    --genomeDir $genomeDir \
    --runThreadN $runThreadN \
    --readFilesIn $readFilesIn_1 $readFilesIn_2 \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattrRGline $outSAMattrRGline \
    --outBAMsortingThreadN $outBAMsortingThreadN \
    --twopassMode "Basic" \
    --outFileNamePrefix $outFileNamePrefix


module purge
module load SAMtools
module list

samtools index -@ $runThreadN $outFileNamePrefix*.bam
