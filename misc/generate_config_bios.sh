#!/bin/bash
# A simple script to generate config file for SbatchAseQuantPipeline

pjdir=~/Documents/projects/wp_ase_dlp
fqdir=$pjdir/inputs/BIOS_RNAseq # prm02 is not mounted on gearshift UI yet, then rsync is an option
idmapfile=$pjdir/inputs/idmapping/BIOS-genotype-rnaseq-ids_with-LL_usable_20210105.txt

mkdir -p $pjdir/scripts/configs/bios

while read -a dnagen_to_rnaseq; do
  if [[ ${dnagen_to_rnaseq[0]} =~ ^'#' ]]; then continue; fi
  uniqId=${dnagen_to_rnaseq[0]}
  sampleId=${dnagen_to_rnaseq[1]}
  fastqId=${dnagen_to_rnaseq[2]}
  cohortId=${dnagen_to_rnaseq[3]}

  metarunfile=$pjdir/scripts/configs/$cohortId-metarun.sh
  if [[ ! -e $metarunfile ]]; then
    echo -e '#Config meta files for '$cohortId > $metarunfile
    echo -e 'set -Eeu -o pipefail' >> $metarunfile
  fi

  conffile=$pjdir/scripts/configs/bios/$uniqId.$fastqId.$sampleId.conf
  cat > "$conffile" <<EOF
# Configuraton file for SbatchWaspPipeline working on $uniqId, $fastqId, $sampleId
fastqId=$fastqId
cohortId=$cohortId
sampleId=$sampleId

# Master directory
projDir=/groups/umcg-bios/tmp01/users/umcg-zzhang/projects/wp_ase_dlp
workDir=\$projDir/outputs/aseQuan_v3

# SLURM logs
logDir=\$workDir/\$cohortId/logDir 

# Temporary directory and files
tmpDir=\$workDir/\$cohortId/tmpDir
fastpTmpDir=\$tmpDir/\$fastqId/fastpTmpDir
starTmpDir=\$tmpDir/\$fastqId/starTmpDir
waspTmpDir=\$tmpDir/\$fastqId/waspTmpDir
byaseTmpDir=\$tmpDir/\$fastqId/byaseTmpDir
gatkTmpDir=\$tmpDir/\$fastqId/gatkTmpDir

# The final output results
optDir=\$workDir/\$cohortId/optDir
fastpOptDir=\$optDir/\$fastqId/fastpOptDir
starOptDir=\$optDir/\$fastqId/starOptDir
waspOptDir=\$optDir/\$fastqId/waspOptDir
byaseOptDir=\$optDir/\$fastqId/byaseOptDir
gatkOptDir=\$optDir/\$fastqId/gatkOptDir

# Genome sequence and gene structure
# NOTE: genomeAnnotationFile is used by STAR (for genome index) and Byase.
genomeDir=\$workDir/genomeDir
genomeFastaFile=\$projDir/inputs/GRCh37_reference/human_g1k_v37.fasta
genomeAnnotationFile=\$projDir/inputs/Gencode/gencode.annotation.ensembl_canonical.protein_coding.GRCh37.gff3.gz

# Genotypes (GT field is required)
snpH5dbDir=\$workDir/snpH5dbDir/\$cohortId
vcfFile=\$projDir/outputs/phasing/all-\$cohortId-singleAlt.vcf.gz
chromInfoFile=\$projDir/inputs/Ensembl_references/human_g1k_v37_chrom_info.txt

# FASTQ files
fastqDir=\$workDir/\$cohortId/tmpDir
fastqPrefix=
fastqSuffix=.fq.gz

# Gene ids TODO: could be removed safely
# geneIdFile=\$projDir/inputs/Ensembl_references/protein_coding_gene_id.txt

# Tools versions, Python virtual env
PythonEnv=~/Documents/projects/wp_ase_dlp/scripts/.env
GCCVer=GCC/7.3.0-2.30
HDF5Ver=HDF5/1.8.14-foss-2018b
STARVer=STAR/2.6.1c-foss-2018b
PythonVer=Python/3.7.4-GCCcore-7.3.0-bare
BCFtoolsVer=BCFtools/1.11-GCCcore-7.3.0
SAMtoolsVer=SAMtools/1.9-foss-2018b
GATKVer=GATK/4.1.4.1-Java-8-LTS  # Only works for GATK/ASEReadCounter
# GSLver=GSL/2.5-GCCcore-7.3.0     # Only works for aScan
# BoostVer=Boost/1.67.0-foss-2018b # Only works for aScan

# vim: set nowrap ft=sh ts=4 tw=120:
EOF
# For BIOS samples you have to rsync files from the old machine (calculon.hpc.rug.nl) to the new one (gearshift.hpc.rug.nl)
  src_dir="umcg-zzhang@172.23.34.247:$fqdir/${fastqId}_*"
  dst_dir="$pjdir/outputs/aseQuan_v3/$cohortId/tmpDir/"
  rsync_cmd="rsync -azqhL $src_dir $dst_dir"
  sbjob_cmd="$pjdir/scripts/bin/SbatchAseQuantPipeline -c $conffile"
  echo "echo -e \$LINENO: && ${rsync_cmd} && ${sbjob_cmd}" >> $metarunfile
done < ${idmapfile}

