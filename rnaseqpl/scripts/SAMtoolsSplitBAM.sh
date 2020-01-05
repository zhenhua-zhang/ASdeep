#!/bin/bash
#
## Split the first STAR mapped reads by chr
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

opt=$(getopt -l "workDir:,fastqId:,help" -- "w:i:h" $@)
eval set -- ${opt}
while true; do
	case $1 in
		-w|--workDir) shift && workDir=$1 ;;
		-i|--fastqId) shift && fastqId=$1 ;;
		--help) echo_help && exit 0 ;;
		--) shift && break;;
	esac
	shift
done

module purge
module load SAMtools
module list

threads=$(check_cpus)
chromId=$SLURM_ARRAY_TASK_ID

waspTmpdir=$workDir/tmpdir/$fastqId/waspTmpdir/$chromId/map1

# Long time sleeping help to avoid conflicts. ;D
sleep $[ $[ 22 - $SLURM_ARRAY_TASK_ID ] * 2 ]

samtools view -hb \
    -@ $threads \
    -o $waspTmpdir/${fastqId}_${chromId}.bam \
    $workDir/tmpdir/$fastqId/starTmpdir/$fastqId*.bam \
    $chromId

samtools sort \
    -@ $threads \
    -o $waspTmpdir/${fastqId}_${chromId}_sort.bam \
    $waspTmpdir/${fastqId}_${chromId}.bam

rm -fr $waspTmpdir/${fastqId}_${chromId}.bam

samtools index \
    -@ $threads \
    $waspTmpdir/${fastqId}_${chromId}_sort.bam

