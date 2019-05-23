#!/bin/bash

module load hisat2

while getopt "r:v:h" OPT; do
	case ${OPT} in 
		r)
			reference=${OPTARG}
			;;
		v)
			variants=${OPTARG}
			;;
		h)
			usage() && exit -1
		\?)
			usage() && exit -1
			;;
	esac
done
shift $((OPTIND-1))

REFERENCE=
VARIANTS=
hisat2_extract_snps_haplotypes_VCF.py
