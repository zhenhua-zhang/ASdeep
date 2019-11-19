#!/bin/bash
#
## Split the first STAR mapped reads by chr
#
set -o errexit
set -o errtrace

opt=$(getopt -l "workDir:,fastqId:,map1Bam:,help" -- "w:s:m:h" $@)
eval set -- ${opt}
while true; do
	case $1 in
		-w|--workDir)
			shift && workDir=$1 ;;
		-s|--fastqId)
			shift && fastqId=$1 ;;
		-m|--map1Bam)
			shift && map1Bam=$1 ;;
		-h|--help) 
			echo_help && exit 0 ;;
		--)
			shift && break;;
	esac
	shift
done

module purge
module load SAMtools
module list

threads=$SLURM_CPUS_PER_TASK
chromId=$SLURM_ARRAY_TASK_ID
map1PerChrDir=$workDir/$fastqId/waspTmpdir/$SLURM_ARRAY_TASK_ID/map1dir

mkdir -p $map1PerChrDir

samtools view -hb -@ $threads -o $map1PerChrDir/$chromId.map1.bam $map1Bam $chromId
samtools sort -@ $threads -o $map1PerChrDir/$chromId.map1.sort.bam $map1PerChrDir/$chromId.map1.bam
samtools index -@ $threads $map1PerChrDir/$chromId.map1.sort.bam

# clean up
# rm $map1PerChrDir/$chromId.map1.bam -fr


