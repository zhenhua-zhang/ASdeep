#!/bin/bash
#
## Filter remapped reads
#
set -o errexit
set -o errtrace

source /apps/modules/modules.bashrc
source ${WASPPL_SCRIPT_PATH}/utils.sh

echo_help() {
    cat <<EOF

Usage:
  bash $(basename $0) [options]

Help:
  -w|--workDir [Required]
      Working directory.
  -i|--fastqId [Required]
      FASTQ files ID.
  --help
      Print this help context and exit.

More information please contact Zhenhua Zhang <zhenhua.zhang217@gmail.com>

EOF
    exit 0
}

long_opts="workDir:,fastqId:,waspPath:,virtualEnv:,chromInfo:,snph5db:,help"
opt=$(getopt -l $long_opts -- "w:i:W:v:c:s:d:" $@)
eval set -- ${opt}
while true; do
	case $1 in
        -w|--workDir) shift && workDir=$1 ;;
        -i|--fastqId) shift && fastqId=$1 ;;
        -W|--waspPath) shift && waspPath=$1 ;;
		-v|--virtualEnv) shift && virtualEnv=$1 ;;
        -c|--chromInfo) shift && chromInfo=$1 ;;
        -d|--snph5db) shift && snph5db=$1 ;;
        -s|--sampleIdFile) shift && sampleIdFile=$1 ;;
		--help) echo_help ;;
		--) shift && break ;;
	esac
	shift
done

workDir=${workDir:?-w/--workDir is required!}
workDir=$(readlink -f $workDir)

waspPath=${waspPath:=~/tools/WASP} # FIXME: when using sbatch, it doesn't work!!
waspPath=$(readlink -f $waspPath)

fastqId=${fastqId:?-i/--fastqId is required!}

snph5db=${snph5db:?-d/--snph5db is required!}
snph5db=$(readlink -f $snph5db || ERRO "Failed to get the path of $snph5db")

chromId=${SLURM_ARRAY_TASK_ID:=1}
perChromDir=$workDir/tmpdir/$fastqId/waspTmpdir/perChrom

#
## Split a huge BAM file by chromosome.
#
module purge
module load SAMtools
module list

samtools view -hb \
    -o $perChromDir/$chromId/${fastqId}_${chromId}.remap.bam \
    $workDir/tmpdir/$fastqId/waspTmpdir/starMapping/*.bam \
    $chromId

#
## Find intersecting SNPs
#
module purge
module load HDF5/1.8.14-foss-2015b Python/3.6.3-foss-2015b
module list

# Load virtual environment if needed, otherwise will use the loaed or system default Python
# interpreter
if [ $virtualEnv"xxx" != "xxx" ]; then
	source $virtualEnv/bin/activate
fi

#
## Conquer
#
toRemapBAM=$perChromDir/$chromId/${fastqId}_${chromId}.to.remap.bam
remapBAM=$perChromDir/$chromId/${fastqId}_${chromId}.remap.bam

filterRemappedReadsPy=$waspPath/mapping/filter_remapped_reads.py
remapKeptBam=$perChromDir/$chromId/${fastqId}_${chromId}.remap.keep.bam

python $filterRemappedReadsPy \
    $toRemapBAM \
    $remapBAM \
    $remapKeptBam

#
## Merge original and remapped BAM
#
module purge
module load SAMtools
module list

keptBAM=$perChromDir/$chromId/${fastqId}_$chromId.keep.bam
keptMergedBAM=$perChromDir/$chromId/${fastqId}_$chromId.keep.merged.bam
keptMergedSortedBAM=${keptMergedBAM/.bam/.sort.bam}

samtoolsThreads=$(check_cpus 1)

samtools merge -f $keptMergedBAM $remapKeptBam $keptBAM
samtools sort -@ $samtoolsThreads -o $keptMergedSortedBAM $keptMergedBAM
samtools index -@ $samtoolsThreads $keptMergedSortedBAM

#
## Remove duplicates
#
rmdupPePy=$waspPath/mapping/rmdup_pe.py
keptMergedSortedRmdupBam=${keptMergedSortedBAM/.bam/.rmdup.bam}

module purge
module load HDF5/1.8.14-foss-2015b Python/3.6.3-foss-2015b
module list

if [ $virtualEnv"xxx" != "xxx" ]; then
	source $virtualEnv/bin/activate
fi

python $rmdupPePy $keptMergedSortedBAM $keptMergedSortedRmdupBam

#
## Sort and index BAM
#
module purge
module load SAMtools
module list

samtoolsThreads=$(check_cpus 1)
keptMergedSortedRmdupSortBam=${keptMergedSortedRmdupBam/.bam/.sort.bam}
samtools sort -@ $samtoolsThreads -o $keptMergedSortedRmdupSortBam $keptMergedSortedRmdupBam
samtools index -@ $samtoolsThreads $keptMergedSortedRmdupSortBam

#
## Get Allele-specific read counts
#
module purge
module load HDF5/1.8.14-foss-2015b Python/3.6.3-foss-2015b
module list

if [ $virtualEnv"xxx" != "xxx" ]; then
	source $virtualEnv/bin/activate
fi

bam2h5py=$waspPath/CHT/bam2h5.py
individual=$(grep -w $fastqId $sampleIdFile | cut -f2)
refAlleleCountsFile=$perChromDir/$chromId/${fastqId}_${chromId}.refAlleleCounts.h5
altAlleleCountsFile=$perChromDir/$chromId/${fastqId}_${chromId}.altAlleleCounts.h5
otherAlleleCountsFile=$perChromDir/$chromId/${fastqId}_${chromId}.otherAlleleCounts.h5
readsCounts=$perChromDir/$chromId/${fastqId}_${chromId}.allCounts.h5
# readsCountsInText=$perChromDir/$chromId/${fastqId}_${chromId}.readCountsInText.txt

python $bam2h5py \
    --chrom $chromInfo \
    --test_chrom $chromId \
    --snp_index $snph5db/$chromId/snps_index.h5 \
    --snp_tab $snph5db/$chromId/snps_tab.h5 \
    --haplotype $snph5db/$chromId/haplotype.h5 \
    --individual $individual \
    --ref_as_counts $refAlleleCountsFile \
    --alt_as_counts $altAlleleCountsFile \
    --other_as_counts $otherAlleleCountsFile \
    --read_counts $readsCounts \
    $keptMergedSortedRmdupSortBam
    # --txt_counts $readsCountsInText \

