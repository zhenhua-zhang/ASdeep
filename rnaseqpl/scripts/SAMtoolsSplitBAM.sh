#!/bin/bash
#
## Split the first STAR mapped reads by chr
#
set -o errexit
set -o errtrace

source /apps/modules/modules.bashrc
source ${WASPPL_SCRIPT_PATH}/utils.sh

long_opts="workDir:,fastqId:,help"
opt=$(getopt -l $long_opts -- "w:i:h" $@)
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

