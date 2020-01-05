#!/bin/bash
#
## Find intersected SNPs
#

# XXX: 
#    1. Using "job array" by Slurm

set -o errexit
set -o errtrace

source /apps/modules/modules.bashrc

# [ $SCRIPT_PATH"xxx" == "xxx" ] \
#     && echo "No SCRIPT_PATH is set." && exit 1 \
#     || source $SCRIPT_PATH/utils.sh

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

Help:
  -w|--workDir [Required]
      The working directory.
  -i|--fastqId [Required]
      The id of fastq files.
  -d|--snph5db [Required]
      The HDF5 database for the WASP.
  -s|--sampleIdFile [Required]
      The sample id in the HDF5 database.
  -W|--waspPath [Optional]
      The path to the WASP.
  -v|--virtualEnv [Optional]
      The python virtual environment for the WASP
  -h|--help
    Print this help context and exit.

More information please contact Zhenhua Zhang <zhenhua.zhang217@gmail.com>

EOF
	exit 0
}

opt=$(getopt -l "workDir:,fastqId:,snph5db:,waspPath:,virtualEnv:,sampleIdFile:,help" -- "w:i:d:W:v:s:h" $@)
eval set -- ${opt}
while true; do
	case $1 in
        -w|--workDir) shift && workDir=$1 ;;
        -i|--fastqId) shift && fastqId=$1 ;;
		-d|--snph5db) shift && snph5db=$1 ;;
        -s|--sampleIdFile) shift && sampleIdFile=$1 ;;
        -W|--waspPath) shift && waspPath=$1 ;;
		-v|--virtualEnv) shift && virtualEnv=$1 ;;
		--help) echo_help ;;
		--) shift && break ;;
	esac
	shift
done

workDir=${workDir:?-w/--workDir is required!}
workDir=$(readlink -f $workDir)

fastqId=${fastqId:?-i/--fastqId is required!}
snph5db=${snph5db:?-d/--snph5db is required!}
sampleIdFile=${sampleIdFile:?-s/--sampleIdFile is required!}

waspPath=${waspPath:=~/tools/WASP} # FIXME: when using sbatch, it doesn't work!!
waspPath=$(readlink -f $waspPath)

chromId=${SLURM_ARRAY_TASK_ID:=1}

#
## Split a huge BAM file by chromosome.
#
module purge
module load SAMtools
module list

perChromDir=$workDir/tmpdir/$fastqId/waspTmpdir/perChrom
mkdir -p $perChromDir/$chromId

samtools view -hb \
    -o $perChromDir/$chromId/${fastqId}_${chromId}.bam \
    $workDir/tmpdir/$fastqId/starTmpdir/*.bam \
    $chromId

samtools index $perChromDir/$chromId/${fastqId}_${chromId}.bam


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

# Load HDF5 database for VCF or imputation results. The database should be constructed before.
snph5db=${snph5db:=$workDir/snph5db}
snpTabFile=$snph5db/snps_tab.h5
snpIndexFile=$snph5db/snps_index.h5
haplotypeFile=$snph5db/haplotype.h5

# find_intersecting_snps.py
findIntersectingSnpsPy=$waspPath/mapping/find_intersecting_snps.py
python $findIntersectingSnpsPy \
    --is_paired_end \
    --is_sorted \
    --output_dir $perChromDir/$chromId \
    --snp_tab $snpTabFile \
    --snp_index $snpIndexFile \
    --haplotype $haplotypeFile \
    --samples $(grep $fastqId $sampleIdFile | cut -f1) \
    $perChromDir/$chromId/${fastqId}_${chromId}.bam

