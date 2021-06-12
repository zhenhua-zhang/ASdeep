#!/bin/bash
# A simple script to generate config file for SbatchAseQuantPipeline

pjdir=~/Documents/projects/wp_ase_dlp
fqdir=${fqdir:=/groups/umcg-bios/prm02/rawdata/rnaseq} # prm02 is not mounted on gearshift UI yet, then rsync is an option
idmapfile=${idmapfile:=$pjdir/inputs/idmapping/BIOS-genotype-rnaseq-ids_with-LL_usable_20210105.txt}

mkdir -p $pjdir/scripts/configs/bios

while read -a dnagen_to_rnaseq; do
    if [[ ${dnagen_to_rnaseq[0]} =~ ^'#' ]]; then continue; fi
    uniqId=${dnagen_to_rnaseq[0]}
    sampleId=${dnagen_to_rnaseq[1]}
    fastqId=${dnagen_to_rnaseq[2]}
    cohortId=${dnagen_to_rnaseq[3]}

    metarunfile=$pjdir/scripts/configs/BIOS-metarun-aseq.sh
    if [[ ! -e $metarunfile ]]; then
        echo -e '#Config meta files for BIOS' > $metarunfile
        echo -e 'set -Eeu -o pipefail' >> $metarunfile
    fi

    conffile=$pjdir/scripts/configs/bios/$uniqId.$fastqId.$sampleId.conf
    cat > "$conffile" <<EOF
# Configuraton file for SbatchWaspPipeline working on $uniqId, $fastqId, $sampleId
fastqId=$fastqId
cohortId=$cohortId
sampleId=$sampleId

# Mater directory
projDir=/groups/umcg-bios/tmp01/users/umcg-zzhang/projects/wp_ase_dlp
workDir=\$projDir/outputs/aseQuan_v2

# SLURM logs
logDir=\$workDir/\$cohortId/logDir 

# Temporary directory and files
tmpDir=\$workDir/\$cohortId/tmpDir
fastpTmpDir=\$tmpDir/\$fastqId/fastpTmpDir
starTmpDir=\$tmpDir/\$fastqId/starTmpDir
waspTmpDir=\$tmpDir/\$fastqId/waspTmpDir
gatkTmpDir=\$tmpDir/\$fastqId/gatkTmpDir
aseqTmpDir=\$tmpDir/\$fastqId/aseqTmpDir

# The final output results
optDir=\$workDir/\$cohortId/optDir
fastpOptDir=\$optDir/\$fastqId/fastpOptDir
starOptDir=\$optDir/\$fastqId/starOptDir
waspOptDir=\$optDir/\$fastqId/waspOptDir
gatkOptDir=\$optDir/\$fastqId/gatkOptDir
aseqOptDir=\$optDir/\$fastqId/aseqOptDir

# Genome sequence and gene structure
genomeDir=\$workDir/genomeDir
genomeFastaFile=\$projDir/inputs/GRCh37_reference/human_g1k_v37.fasta
genomeAnnotationFile=\$projDir/inputs/Ensembl_references/Homo_sapiens.GRCh37.75.gtf

# Genotypes (GT field is required)
snpH5dbDir=\$workDir/snpH5dbDir/\$cohortId
vcfFile=\$projDir/outputs/phasing/all-\$cohortId-singleAlt.vcf.gz
chromInfoFile=\$projDir/inputs/Ensembl_references/human_g1k_v37_chrom_info.txt

# FASTQ files
fastqDir=\$workDir/\$cohortId/tmpDir
fastqPrefix=
fastqSuffix=.fq.gz

# Gene ids
geneIdFile=\$projDir/outputs/aseReports/summary/BIOS_ASE_genes.txt

# Tools versions, Python virtual env
pyEnv=~/Documents/projects/wp_ase_dlp/scripts/.env
GCCVer=GCC/7.3.0-2.30
GATKVer=GATK/4.1.4.1-Java-8-LTS
HDF5Ver=HDF5/1.8.14-foss-2018b
STARVer=STAR/2.6.1c-foss-2018b
PythonVer=Python/3.7.4-GCCcore-7.3.0-bare
BCFtoolsVer=BCFtools/1.11-GCCcore-7.3.0
SAMtoolsVer=SAMtools/1.9-foss-2018b

# vim: set nowrap ft=sh ts=4 tw=120:
EOF
    # For BIOS samples you have to rsync files from the old machine (calculon.hpc.rug.nl) to the
    # new one (gearshift.hpc.rug.nl)
    # echo "echo -ne \$LINENO: && rsync -avzh umcg-zzhang@172.23.34.247:$fqdir/${fastqId}_* $pjdir/outputs/aseQuan_v2/$cohortId/tmpDir/ && $pjdir/scripts/bin/SbatchAseQuantPipeline -c $conffile" >> $metarunfile
    echo "echo -ne \$LINENO: && $pjdir/scripts/bin/SbatchAseQuantPipeline-aseq -c $conffile" \
        >> $metarunfile
done < ${idmapfile}
