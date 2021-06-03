#!/bin/bash
# A simple script to generate config file for SbatchAseQuantPipeline

projDir=~/Documents/projects/wp_ase_dlp
idmapfile=$projDir/inputs/geuvadis/samples/igsr_geuvadis-genotype-rnaseq-ids_20210308.txt
configDir=$projDir/scripts/configs/geuvadis

baseUrl=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR188

mkdir -p $configDir
while read -a myarray; do
    if [[ ${myarray[0]} =~ ^'#' ]]; then continue; fi
    uniqId=${myarray[0]}
    sampleId=${myarray[1]}
    fastqId=${myarray[2]}
    cohortId=${myarray[3]}

    metarunfile=$projDir/scripts/configs/$cohortId-metarun.sh
    if [[ ! -e $metarunfile ]]; then
        echo -e '#Config meta files for '$cohortId > $metarunfile
        echo -e 'set -Eeu -o pipefail' >> $metarunfile
    fi

    conffile=$configDir/$uniqId.$fastqId.$sampleId.conf
    cat > "$conffile" <<EOF
# Configuraton file for SbatchWaspPipeline working on $uniqId, $fastqId, $sampleId
fastqId=$fastqId
cohortId=$cohortId
sampleId=$sampleId

# Mater directory
projDir=$projDir
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
vcfFile=\$projDir/outputs/phasing/ALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.deploid.snp.nondup.vcf.gz
chromInfoFile=\$projDir/inputs/Ensembl_references/human_g1k_v37_chrom_info.txt

# FASTQ files
fastqDir=\$workDir/\$cohortId/tmpDir
fastqPrefix=
fastqSuffix=.fq.gz

# Gene ids
geneIdFile=\$projDir/inputs/Ensembl_references/protein_coding_gene_id.txt

# Tools versions, Python virtual env
pyEnv=~/Documents/projects/wp_ase_dlp/scripts/.env
GCCVer=GCC/7.3.0-2.30
GATKVer=GATK/4.1.0.0-Java-1.8
HDF5Ver=HDF5/1.10.2-intel-2018b
STARVer=STAR/2.6.0c-foss-2018a
PythonVer=Python/3.7.4-GCCcore-8.3.0
BCFtoolsVer=BCFtools/1.10.2-GCC-9.3.0
SAMtoolsVer=SAMtools/1.10-GCC-9.3.0

# vim: set nowrap ft=sh ts=4 tw=120:
EOF
    
    # For Geuvadis samples, ought to download their FASTQ files from ftp://ftp.sra.ebi.ac.uk/vol1/fastq
    echo "echo -ne \$LINENO: " \
         "&& curl -C- $baseUrl/$fastqId/${fastqId}_1.fastq.gz" \
         "--output $projDir/outputs/aseQuan_v2/$cohortId/tmpDir/${fastqId}_R1.fq.gz" \
         "&& curl -C- $baseUrl/$fastqId/${fastqId}_2.fastq.gz" \
         "--output $projDir/outputs/aseQuan_v2/$cohortId/tmpDir/${fastqId}_R2.fq.gz" \
         "&& $projDir/scripts/bin/SbatchAseQuantPipeline -c $conffile" >> $metarunfile
done < ${idmapfile}
