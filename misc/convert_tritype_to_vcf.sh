#!/bin/bash
awk 'BEGIN{count=0} {if ($1 == LAST_LINE) {count++; print $1"."count} else {count=0; LAST_LINE=$1; print $1}}' SNPs.txt > SNPs_mkdup.txt
awk 'BEGIN{count=0} {if ($3 == LAST_LINE) {count++; print $0"."count} else {count=0; LAST_LINE=$3}}' SNPMappings.txt > SNPMappings_mkdup.txt

module load GenotypeHarmonizer

java -jar $EBROOTGENOTYPEHARMONIZER/GenotypeHarmonizer.jar \
    -i ${trityper_dir} \
    -I TRITYPER \
    -o ${trityper_dir##-TriTyper} \
    -O PLINK_BED

