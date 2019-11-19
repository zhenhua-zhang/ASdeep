#!/bin/bash
#
## STAR alignReads. First round mapping FASTQ reads to reference genome.
#

# TODO:
# 	1. A better allocation of threads for runThreadN and outBAMsortingThreadN
# 	2. Specify genomeDir

module purge
module load STAR/2.5.1b-foss-2015b
module list

runThreadN=$[ $SLURM_CPUS_PER_TASK - 2 ]
outBAMsortingThreadN=2

readFilesIn_1=$fastq_paired_1
readFilesIn_2=$fastq_paired_2
readFilesCommand="zcat"  # If the fataq or fasta files are compressed the --readFilesCommand options

outSAMtype_word_1=BAM  # The first word for outSAMtype
outSAMtype_word_2=SortedByCoordinate  # The second word for outSAMtype

outSAMattrRGline="ID:$fastq_id SM:$fastq_id PL:ILLUMINA"

twopassMode="Basic"  # Two-pass mode

map1_outFileNamePrefix=$star_tmdir/$fastq_id

STAR \
    --runMode alignReads \
    --genomeDir $genomeDir \
    --runThreadN $runThreadN \
    --readFilesIn $readFilesIn_1 $readFilesIn_2 \
	--readFilesCommand $readFilesCommand \
    --outSAMtype $outSAMtype_word_1 $outSAMtype_word_2 \
    --outSAMattrRGline $outSAMattrRGline \
    --outBAMsortingThreadN $outBAMSortingThreadN \
    --twopassMode $twopassMode \
    --outFileNamePrefix $map1_outFileNamePrefix


module purge
module load SAMtools
module list

samtools index $map1_bam
