#!/bin/bash
#
## STAR index. Construction of the genome index for STAR aligner of RNA-seq data. 
#

source ../bin/utils

echo_help() {
	cat <<EOF

Help:
  -w|--workDir Required
      Working directory
  -g|--genomeDir Required
  -f|--genomeFastaFile Required
  -a|--genomeAnnotationsFile Required
  --runThreadN Optional
      Number of threadings STAR to build the genome index
  --help Optional
      Print this help context and exit.

More information please contact Zhenhua Zhang <zhenhua.zhang217@gmail.com>
EOF
	exit 0
}

opt=$(getopt -l "workDir:,genomeDir:,genomeFastaFile:,genomeAnnotationsFile:,runThreadN:,help" -- "w:f:g:a:" $@)
eval set -- ${opt}
while true; do
	case $1 in
		-w|--workDir) shift && workdir=$1 ;;
		-g|--genomeDir) shift && genomeDir=$1 ;;
		-f|--genomeFastaFile) shift && genomeFastaFile=$1;;
		-a|--genomeAnnotationsFile) shift && genomeAnnotationsFile=$1 ;;
		--runThreadN) shift && runThreadN=$1 ;;
		--outFileNamePrefix) shift && outFileNamePrefix=$1 ;;
		--help) echo_help ;;
		--) echo_help;;
	esac
	shift
done

# Variables by command line
workDir=${workDir:?-w/--workDir is required!}
genomeDir=${workDir}/${genomeDir:?-g/--genomeDir is required}
genomeFastaFiles=${genomeFastaFiles:=/apps/data/ftp.broadinstitute.org/bundle/2.8/b37/human_g1k_v37.fasta}
genomeAnnotationsFile=${genomeAnnotationsFile:=/apps/data/ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf}
runThreadN=${runThreadN:=$SLURM_CPUS_PER_TASK}  # Threadings. #FIXME: Only using SLURM
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
