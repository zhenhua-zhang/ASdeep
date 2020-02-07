#!/bin/bash
# Collect outputs in the pipeline to the final output directory

set -o errexit
set -o errtrace

source /apps/modules/modules.bashrc

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

Help:
  -w|--workDir [Required]
      Working directory.
  -i|--fastqId [Required]
      FASTQ files ID.
  -h|--help
      Print this help context and exit.

More information please contact Zhenhua Zhang <zhenhua.zhang217@gmail.com>

EOF
    exit 0
}

opt=$(getopt -l "workDir:,fastqId:,help" -- "w:i:h" $@)
eval set -- $opt

while true; do
    case $1 in
        -w|--workDir) shift && workDir=$1 ;;
        -i|--fastqId) shift && fastqId=$1 ;;
        -h|--help) shift && echo_help ;;
        --) shift && break ;;
    esac

    shift
done

workDir=${workDir:?[E]: -w/--workDir is required!}
workDir=$(readlink -f $workDir)

fastqId=${fastqId:?[E]: -i/--fastqId is required!}

mkdir -p $workDir/optdir/$fastqId/{fastp,star,wasp}Optdir

cp -fr $workDir/tmpdir/$fastqId/fastpTmpdir/*report{.html,.json} $workDir/optdir/$fastqId/fastpOptdir
cp -fr $workDir/tmpdir/$fastqId/starTmpdir/* $workDir/optdir/$fastqId/starOptdir

threads=$(check_cpus)
~/tools/bin/parallel -j $threads \
    mkdir $workDir/optdir/$fastqId/waspOptdir/perChrom/{1} \;
    cp -fr $workDir/tmpdir/$fastqId/waspTmpdir/perChrom/{1}/*.h5 \
    $workDir/optdir/$fastqId/waspOptdir/perChrom/{1} \
    ::: {1..22}

# cat $workDir/tmpdir/$fastqId/waspTmpdir/perChrom/*/*.readCountsInText.txt \
#     > $workDir/optdir/$fastqId/waspOptdir/$fastqId.readCountsInText.txt

module purge
module load SAMtools
module list

samtools merge -fcp \
    -@ $threads \
    $workDir/optdir/$fastqId/waspOptdir/$fastqId.keep.merge.rmdup.bam \
    $workDir/tmpdir/$fastqId/waspTmpdir/perChrom/*/*.keep.merged.sort.rmdup.sort.bam

samtools sort \
    -@ $threads \
    -o $workDir/optdir/$fastqId/waspOptdir/$fastqId.keep.merge.rmdup.sort.bam \
    $workDir/optdir/$fastqId/waspOptdir/$fastqId.keep.merge.rmdup.bam

samtools index \
    -@ $threads \
    $workDir/optdir/$fastqId/waspOptdir/$fastqId.keep.merge.rmdup.sort.bam

rm -fr $workDir/optdir/$fastqId/waspOptdir/$fastqId.keep.merge.rmdup.bam
rm -fr $workDir/tmpdir/$fastqId
