#!/bin/bash
# A simple script to generate config file for SbatchAseQuantPipeline

wkdir="$HOME/Documents/projects/wp_ase_dlp"
idmap=$wkdir/inputs/Geuvadis/samples/igsr_geuvadis-genotype-rnaseq-ids_20210308.txt
confdir=$wkdir/scripts/configs/Geuvadis

baseUrl=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR188

mkdir -p $confdir
while read -a myarray; do
    if [[ ${myarray[0]} =~ ^'#' ]]; then continue; fi
    uniqId=${myarray[0]}
    sampleId=${myarray[1]}
    fastqId=${myarray[2]}
    cohortId=${myarray[3]}

    metarunfile=$wkdir/scripts/configs/$cohortId-metarun.sh
    if [[ ! -e $metarunfile ]]; then
        echo -e '#Config meta files for '$cohortId > $metarunfile
        echo -e 'set -Eeu -o pipefail' >> $metarunfile
    fi

    conffile=$confdir/$uniqId.$fastqId.$sampleId.conf
    cat > "$conffile" <<EOF
# Configuraton file for SbatchWaspPipeline working on $uniqId, $fastqId, $sampleId
fastqId=$fastqId
cohortId=$cohortId
sampleId=$sampleId

# Master directory
projDir=/groups/umcg-bios/tmp01/users/umcg-zzhang/projects/wp_ase_dlp
workDir=\$projDir/outputs/allelic_read_counts

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
vcfFile=\$projDir/outputs/haplotypes/all-\$cohortId-singleAlt.vcf.gz
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
    
    # For Geuvadis samples, ought to download their FASTQ files from ftp://ftp.sra.ebi.ac.uk/vol1/fastq
    echo "echo -ne \$LINENO:" \
         "&& curl -C- --fail"\
         "-o '$wkdir/outputs/allelic_read_counts/$cohortId/tmpDir/${fastqId}_R#1.fq.gz'" \
         "'$baseUrl/$fastqId/${fastqId}_[1-2].fastq.gz'" \
         "&& $wkdir/scripts/misc/SbatchAseQuantPipeline -c $conffile" >> $metarunfile
done < ${idmap}
