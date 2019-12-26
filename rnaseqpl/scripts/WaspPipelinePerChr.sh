#!/bin/bash
#
## Find intersected SNPs, which could have mapping biases
#

set -o errexit
set -o errtrace

echo_help() {
	cat <<EOF

Help:
  -h, --help    Optional. Action: print_info
    Print this help context and exit.

More information please contact Zhenhua Zhang <zhenhua.zhang217@gmail.com>
EOF
	exit 0
}

opt=$(getopt -l "workDir:,FASTQID:,chromInfo:,firstMapBAM:,snpHDF5db:,help" -- "w:s:c:f:H:h" $@)
eval set -- ${opt}
while true; do
	case $1 in
		-c|--chromInfo) shift && chromInfo=$1 ;;
		-f|--map1Bam) shift && map1Bam=$1 ;;
        -w|--waspPath) shift && waspPath=$1 ;;
		-H|--snpHDF5db) shift && snpHDF5db=$1 ;;
		-v|--virtualEnv) shift && virtualEnv=$1 ;;
        --outBAMsortingThreadN) shift && outBAMsortingThreadN=$1 ;;
		--help) echo_help ;;
		--) echo -e "The CLI options are empty" && echo_help ;;
	esac
	shift
done

chromId=$SLURM_ARRAY_TASK_ID
waspPerChrDir=$WORKDIR/$FASTQID/waspTmpdir/$chromId
map1Dir=$waspPerChrDir/map1Dir
waspPath=${waspPath:=../bin/WASP}

module purge
module load HDF5/1.8.14-foss-2015b Python/3.6.3-foss-2015b
module list

# Load virtual environment if needed, otherwise will use the loaed or system default Python
# interpreter
if [ -n $virtualEnv ]; then
	source $virtualEnv/bin/activate
fi

# find_intersecting_snps.py
findIntersectingSnpsDir=$waspPerChrDir/findIntersectingSnpsDir
mkdir -p $findIntersectingSnpsDir

# Load HDF5 database for VCF or imputation results. The database should be constructed before.
snpHDF5db=${snpHDF5db:=${workDir}/snpHDF5db}
snpTabFile=$snpHDF5db/snpsTab.h5
snpIndexFile=$snpHDF5db/snpsIndex.h5
haplotypeFile=$snpHDF5db/haplotype.h5

findIntersectingSnpsPy=${waspPath}/mapping/find_intersecting_snps.py
python $findIntersectingSnpsPy \
    --is_paired_end \
    --is_sorted \
    --output_dir $findIntersectingSnpsDir \
    --snp_tab $snpTabFile \
    --snp_index $snpIndexFile \
    --haplotype $haplotypeFile \
    $map1Dir/$chromId.map1.sort.bam

#
## Remap reads with mapping biases
#
map2dir=$waspPerChrDir/map2dir
mkdir -p $map2_dir

map2ReadFilesIn_1= #TODO: Not yet!!!!
map2ReadFilesIn_2= #TODO: Not yet!!!!

outBAMsortingThreadN=${outBAMsortingThreadN:=2}
[ -n ${SLURM_CPUS_PER_TASK} ] \
    && runThreadN=$[ ${SLURM_CPUS_PER_TASK} - ${outBAMsortingThreadN} ]
    || runThreadN=$[ $(grep -c processor /proc/cpuinfo) - ${outBAMsortingThreadN} ]

module purge
module load STAR/2.5.1b-foss-2015b
module list

STAR \
    --runMode alignReads \
    --genomeDir $genomeDir \
    --runThreadN $runThreadN \
    --readFilesIn $map2ReadFilesIn_1 $map2ReadFilesIn_2 \
    --outSAMattrRGline ID:$FASTQID SM:$FASTQID PL:ILLUMINA \
    --outBAMsortingThreadN $outBAMSortingThreadN \
    --outSAMtype BAM SortedByCoordinate \
    --twopassMode Basic \
	--readFilesCommand zcat \
    --outFileNamePrefix $FASTQID

#
## Filter remapped reads
#
filterRemappedReadsPy=$waspPath/mapping/filter_remapped_reads.py
filterRemappedReadsDir=$waspPerChrDir/filterRemappedReadsDir
mkdir -p $filterRemappedReadsDir

toRemapBAM=
remapBAM=
keptBAM1=

module purge
module load HDF5/1.8.14-foss-2015b Python/3.6.3-foss-2015b
module list

python $filterRemappedReadsPy \
    $toRemapBAM \
    $remapBAM \
    $keptBAM1

#
## Merge original and remapped BAM
#
module purge
module load SAMtools
module list

keptBAM2=
keptMergedBAM=
keptMergedSortedBAM=

samtoolsThreads=$SLURM_CPUS_PER_TASK

samtools merge $keptMergedBAM $keptBAM1 $keptBAM2
samtools sort -@ $samtoolsThreads -o $keptMergedSortedBAM $keptMergedBAM
samtools index -@ $samtoolsThreads $keptMergedSortedBAM

#
## Remove duplicates
#
removeDuplicatedPairedDir=
rmdupPePy=rmdup_pe.py
rmdupBAM=

mkdir -p $removeDuplicatedPairedDir
module purge
module load HDF5/1.8.14-foss-2015b Python/3.6.3-foss-2015b
module list

[ -n $virtualEnv ] && source $virtualEnv/bin/activate
python $rmdupPePy $keptMergedSortedBAM $rmdupBAM

module purge
module load SAMtools
module list

rmdupSotedBAM=
samtools sort -@ $samtoolsThreads -o $rmdupSotedBAM $rmdupBAM
samtools index -@ $samtoolsThreads $rmdupSotedBAM


#
## Get Allele-specific read counts
#
aseDir=$waspTmpdir/aseDir
mkdir -p $aseDir

module purge
module load HDF5/1.8.14-foss-2015b
module load Python/3.6.3-foss-2015b
module list

[ -n $virtualEnv ] && source $virtualEnv/bin/activate

bam2H5Py=$waspPath/CHT/bam2h5.py
individualId=
referenceAlleleCountsFile=${fastqId}_${chromId}_refAlleleCounts.h5
alternativeAlleleCountsFile=${fastqId}_${chromId}_altAlleleCounts.h5
otherAlleleCountsFile=${fastqId}_${chromId}_otherAlleleCounts.h5
readsCounts=${fastqId}_${chromId}_readCounts.h5
readsCountsInText=${fastqId}_${chromId}_readCountsInText.txt

python $bam2H5Py \
    --chrom $chromInfo \
    --snp_index $snpIndexFile \
    --snp_tab $snpTabFile \
    --haplotype $haplotypeFile \
    --individual $individualId \
    --ref_as_counts $referenceAlleleCountsFile \
    --alt_as_counts $alternativeAlleleCountsFile \
    --other_as_counts $otherAlleleCountsFile \
    --read_counts $readsCounts \
    --txt_counts $readsCountsInText \
    $rmdupSotedBAM
