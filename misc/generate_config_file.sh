#!/bin/bash

# idmapfile=${idmapfile:=../../../misc/freeze2_GoNL_related_GTE_20200211_QCpassed.csv}


pjdir=/groups/umcg-bios/tmp01/users/umcg-zzhang/projects/wp_ase_dlp
fqdir=/groups/umcg-bios/prm02/rawdata/rnaseq # prm02 is not mounted on gearshift UI yet, then rsync is an option
idmapfile=${idmapfile:=$pjdir/inputs/idmapping/BIOS-genotype-rnaseq-ids_with-LL_usable_20210105.txt}

mkdir -p $pjdir/scripts/configs/bios

while read -r dnagen_to_rnaseq; do
    if [[ ${dnagen_to_rnaseq:0:1} == "#" ]]; then continue; fi
    sample_id=$(cut -f1 -d$'\t' <<<"$dnagen_to_rnaseq")
    dnavar_id=$(cut -f2 -d$'\t' <<<"$dnagen_to_rnaseq")
    rnaseq_id=$(cut -f3 -d$'\t' <<<"$dnagen_to_rnaseq")
    cohort_id=$(cut -f4 -d$'\t' <<<"$dnagen_to_rnaseq")

    metarunfile=$pjdir/scripts/configs/$cohort_id-metarun.sh
    if [[ ! -e $metarunfile ]]; then
        echo -e '#Config meta files for '$cohort_id > $metarunfile
        echo -e 'set -Eeu -o pipefail' >> $metarunfile
    fi

    conffile=$pjdir/scripts/configs/bios/$rnaseq_id.$sample_id.$dnavar_id.conf
    cat > "$conffile" <<EOF
# Configuraton file for SbatchWaspPipeline working on $sample_id
workDir=$pjdir/outputs/aseQuan
fastqId=$rnaseq_id
fastqDir=\$workDir/tmpDir
fastqPrefix=""
fastqSuffix=.fq.gz

# fastp
fastpExe=~/tools/bin/fastp

# STAR
genomeDir=\$workDir/genomeDir
genomeFastaFile=/apps/data/ftp.broadinstitute.org/bundle/2.8/b37/human_g1k_v37.fasta
genomeAnnotationFile=/apps/data/ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf

# WASP::snp2h5
snp2h5Exe=~/tools/bin/snp2h5
snpH5dbPath=\$workDir/snpH5db/$cohort_id
vcfFileDir=$pjdir/outputs/phasing/${cohort_id}/beagle
vcfFilePrefix=beagle-${cohort_id}-chr
vcfFileSuffix=.vcf.gz
chromInfoFile=$pjdir/misc/human_g1k_v37_chrom_info.txt

# WASP::fasta2h5
fasta2h5Exe=~/tools/bin/fasta2h5
fastaH5dbPath=\$workDir/fastaH5db

# WASP::main
waspDir=~/tools/WASP
pyEnv=$pjdir/scripts/.env
sampleIdFile=$idmapfile
EOF
    echo "rsync -avzh umcg-zzhang@172.23.34.247:$fqdir/${rnaseq_id}_* $pjdir/outputs/aseQuan/tmpDir/ && $pjdir/scripts/bin/SbatchWaspPipeline -c $conffile" >> $metarunfile
done < ${idmapfile}
