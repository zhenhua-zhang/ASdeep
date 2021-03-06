#!/bin/bash
pjdir=/groups/umcg-bios/tmp01/users/umcg-zzhang/projects/wp_ase_dlp
fqdir=/groups/umcg-bios/prm02/rawdata/rnaseq # prm02 is not mounted on gearshift UI yet, then rsync is an option
idmapfile=${idmapfile:=$pjdir/inputs/idmapping/BIOS-genotype-rnaseq-ids_with-LL_usable_20210105.txt}

mkdir -p $pjdir/scripts/configs/bios

while read -r dnagen_to_rnaseq; do
    if [[ ${dnagen_to_rnaseq:0:1} == "#" ]]; then continue; fi
    uniqId=$(cut -f1 -d$'\t' <<<"$dnagen_to_rnaseq")
    sampleId=$(cut -f2 -d$'\t' <<<"$dnagen_to_rnaseq")
    fastqId=$(cut -f3 -d$'\t' <<<"$dnagen_to_rnaseq")
    cohortId=$(cut -f4 -d$'\t' <<<"$dnagen_to_rnaseq")

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
geneIdFile=\$projDir/inputs/Ensembl_references/protein_coding_gene_id.txt

# Python virtual environment
pyEnv=~/Documents/projects/wp_ase_dlp/scripts/.env
EOF
    echo "rsync -avzh umcg-zzhang@172.23.34.247:$fqdir/${fastqId}_* $pjdir/outputs/aseQuan_v2/$cohortId/tmpDir/ && $pjdir/scripts/bin/SbatchAseQuantPipeline -c $conffile" >> $metarunfile
done < ${idmapfile}
