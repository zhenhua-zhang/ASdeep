#!/bin/bash

id_to_id_file=${id_to_id_file:=../../../misc/freeze2_GoNL_related_GTE_20200211_QCpassed.csv}

for id_to_id in $(awk '{print $1";"$2}' ${id_to_id_file});do
    fastq_id=${id_to_id%%;*}
    gonl_id=${id_to_id##*;}
    cat > ../configs/${fastq_id}_${gonl_id}_config <<EOF
# Configuraton file for SbatchWaspPipeline working on ${fastq_id}
projdir=/groups/umcg-bios/tmp04/umcg-zzhang/projects/ASEDLP
workDir=\${projdir}/workdir

# fastp
fastqId=${fastq_id}
fastqDir=\${projdir}/inputs/BIOS_FASTQs

# STAR
genomeDir=\$workDir/genomedir
genomeFastaFile=/apps/data/ftp.broadinstitute.org/bundle/2.8/b37/human_g1k_v37.fasta
genomeAnnotationsFile=/apps/data/ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf

# WASP::snp2h5
snph5db=\${workDir}/snph5db
snp2h5Exe=~/tools/bin/snp2h5
vcfFileDir=\$projdir/inputs/GoNl_genotypes
chromInfoFile=\${projdir}/misc/human_g1k_v37_chrom_info.txt

# WASP::fasta2h5
fastah5db=\${workDir}/fastah5db
fasta2h5Exe=~/tools/bin/fasta2h5

# WASP::MAIN
waspPath=~/tools/WASP-0.3.4
virtualEnv=\${WASPPL_SCRIPT_PATH}/../.env
sampleIdFile=\${projdir}/misc/freeze2_GoNL_related_GTE_30092016_QCpassed.csv
EOF

done
