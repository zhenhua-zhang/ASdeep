#!/bin/bash
# A simple script to generate config file for SbatchAseQuantPipeline

tissue=WHLBLD
tissue=MSCLSK
pjdir=~/Documents/projects/wp_ase_dlp
idmapfile=${idmapfile:-$pjdir/inputs/GTEx/misc/idmapfile-$tissue.tsv}
geneIdFile=$pjdir/outputs/aseReports/GTEx/$tissue/asedlp_report_exMHC_candidate_genes.txt

if [[ ! -f $pjdir/inputs/GTEx/misc/idmapfile-$tissue.tsv ]]; then
    echo -e "#uniq_id\tdna_id\trna_id\tbiobank_id" > $pjdir/inputs/GTEx/misc/idmapfile-$tissue.tsv
    for x in $pjdir/inputs/GTEx/ase/*.gz; do
        if [[ $(zgrep -cm 1 $tissue ${x}) -eq 1 ]]; then
            x=$(basename ${x})
            x=${x%%.*}
            echo -e $x"\t"$x"\t"$x"\tGTEx"
        fi
    done >> $pjdir/inputs/GTEx/misc/idmapfile-$tissue.tsv
fi

mkdir -p $pjdir/scripts/configs/GTEx-$tissue
while read -a dnagen_to_rnaseq; do
    if [[ ${dnagen_to_rnaseq[0]} =~ ^'#' ]]; then continue; fi
    uniqId=${dnagen_to_rnaseq[0]}
    sampleId=${dnagen_to_rnaseq[1]}
    fastqId=${dnagen_to_rnaseq[2]}
    cohortId=${dnagen_to_rnaseq[3]}

    metarunfile=$pjdir/scripts/configs/$cohortId-$tissue-metarun.sh
    if [[ ! -e $metarunfile ]]; then
        echo -e '#Config meta files for '$cohortId > $metarunfile
        echo -e 'set -Eeu -o pipefail' >> $metarunfile
    fi

    conffile=$pjdir/scripts/configs/GTEx-$tissue/$uniqId.$fastqId.$sampleId.conf
    cat > $conffile <<EOF
# Configuraton file for SbatchWaspPipeline working on $uniqId, $fastqId, $sampleId
fastqId=$fastqId
cohortId=$cohortId
sampleId=$sampleId

# Mater directory
projDir=~/Documents/projects/wp_ase_dlp
workDir=\$projDir/outputs/allelic_read_counts

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
genomeAnnotationFile=\$projDir/inputs/Ensembl_references/Homo_sapiens.GRCh37.75.db

# Genotypes (GT field is required)
snpH5dbDir=\$workDir/snpH5dbDir/\$cohortId
vcfFile=\$projDir/outputs/haplotypes/gtex/phASER_WASP_GTEx_v8_merged_hq_b37_annot.vcf.gz
chromInfoFile=\$projDir/inputs/Ensembl_references/human_g1k_v37_chrom_info.txt

# FASTQ files
fastqDir=\$workDir/\$cohortId/tmpDir
fastqPrefix=
fastqSuffix=.fq.gz
tissue=$tissue

# Gene ids
geneIdFile=$geneIdFile

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
    echo "echo -ne \$LINENO: && $pjdir/scripts/bin/SbatchAseQuantPipeline-gtex -c $conffile" >> $metarunfile
done < ${idmapfile}

echo "# vim: set nowrap number:" >> $metarunfile
