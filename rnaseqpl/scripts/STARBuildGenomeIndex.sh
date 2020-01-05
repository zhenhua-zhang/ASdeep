#!/bin/bash
#
## STAR index. Construction of the genome index for STAR aligner of RNA-seq data. 
#

# TODO:
#   1. --genomeAnnotationsFile should be optional. but now it's required.

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

Help:
  -w|--workDir [Required]
      Working directory.
  -g|--genomeDir [Required]
      Directory for the indexed genome.
  -f|--genomeFastaFile [Required]
      Fasta file of the genome.
  -a|--genomeAnnotationsFile [Required]
      GFF/GTF annotation file. Should be UNZIPPED!
  --runThreadN [Optional]
      Number of threadings for STAR to build the genome index.
  --outFileNamePrefix [Optional]
      A prefix string for the genome index.
  --help
      Print this help context and exit.

More information please contact Zhenhua Zhang <zhenhua.zhang217@gmail.com>

EOF
	exit 0
}

opt=$(getopt -l "workDir:,genomeDir:,genomeFastaFile:,genomeAnnotationsFile:,runThreadN:,outFileNamePrefix:,help" -- "w:f:g:a:" $@)
eval set -- $opt
while true; do
	case $1 in
		-w|--workDir) shift && workDir=$1 ;;
		-g|--genomeDir) shift && genomeDir=$1 ;;
		-f|--genomeFastaFile) shift && genomeFastaFile=$1;;
		-a|--genomeAnnotationsFile) shift && genomeAnnotationsFile=$1 ;;
		--runThreadN) shift && runThreadN=$1 ;;
		--outFileNamePrefix) shift && outFileNamePrefix=$1 ;;
		--help) echo_help ;;
		--) shift && break ;;
	esac
	shift
done

# Variables by command line
workDir=${workDir:?-w/--workDir is required!}
genomeDir=${genomeDir:?-g/--genomeDir is required}
genomeFastaFiles=${genomeFastaFiles:=/apps/data/ftp.broadinstitute.org/bundle/2.8/b37/human_g1k_v37.fasta}
genomeAnnotationsFile=${genomeAnnotationsFile:=/apps/data/ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf}
runThreadN=${runThreadN:=$(check_cpus)}
outFileNamePrefix=${outFileNamePrefix:=STARGenomeIndex}

module purge
module load STAR/2.5.1b-foss-2015b
module list

STAR \
  --runMode genomeGenerate \
  --genomeDir $genomeDir \
  --genomeFastaFiles $genomeFastaFiles \
  --sjdbGTFfile $genomeAnnotationsFile \
  --runThreadN $runThreadN \
  --outFileNamePrefix $outFileNamePrefix
