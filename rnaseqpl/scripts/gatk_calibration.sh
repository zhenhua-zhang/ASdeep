#!/bin/bash
#SBATCH --time=:10:0
#SBATCH --mem=5G
#SBATCH --cpus=10
#SBATCH --output=%j-%u-gatk_calibration.log
#SBATCH --job-name=gatk_calibration

set -o errexit
set -o errtrace
source /apps/modules/modules.bashrc

# spdir=$(dirname $(readlink -f "$0"))
# source "$spdir"/presettings.sh
source presettings.sh

#
## GATK: Base quality score recalibration
#

module load GATK/4.1.2.0-Java-1.8.0_144-unlimited_JCE

#
## TODO: MergeBamAlignment
#

#
## MarkDuplicates
#
$gatk_exe \
    MarkDuplicates \
    --INPUT $input_bam \
    --OUTPUT $sort_mkdp_bam \
    --CREATE_INDEX true \
    --METRICS_FILE ${sort_mkdp_metrics} \
    --VALIDATIONG_STRINGENCY SILENT


#
## SplitNCigarReads
#
$gatk_exe \
    SPlitNCigarReads \
    -R $rfr_genome \
    -I $sort_mkdp_bam \
    -O $sort_mkdp_sncr_bam


#
## BaseRecalibrator
#
$gatk_exe \
    --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails -Xloggc:gc_log.log -Xms4000m" \
    BaseRecalibrator \
    -R $rfr_genome \
    -I $sort_mkdp_sncr_bam \
    -O $sort_mkdp_sncr_bsqr_rcl \
    -known-sites $rfr_snps \
    -known-sites $rfr_indels \
    --use-original-qualities

$gatk_exe \
    --java-options "-XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails -Xloggc:gc_log.log -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms3000m" \
    ApplyBQSR \
    -R $rfr_genome \
    -I $sort_mkdp_sncr_bam \
    -O $sort_mkdp_sncr_bsqr_bam \
    --bqsr-recal-file $sort_mkdp_sncr_bsqr_rcl \
    --add-output-sam-program-record \
    --use-original-qualities


# #
# ## Create an index file to 
# #
# module load SAMtools

# $smtl_exe index \
#     -@ $smtl_thd \
#     $sort_mkdp_sncr_bsqr_bam \
#     $sort_mkdp_sncr_bsqr_bai
