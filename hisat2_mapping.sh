#!/bin/bash
#SBATCH --time=3:0:0
#SBATCH --mem=15G
#SBATCH --cpus=15
#SBATCH --output=%j-%u-hisat_mapping.log
#SBATCH --job-name=hisat_mapping

DEBUG=0

ERRO() {
	echo -e "[ERRO]: $1" >&2 && exit -1
}

WARN() {
	echo -e "[WARN]: $1" >&2
}

INFO() {
	echo -e "[INFO]: $1"
}

check_dir() {
    for DIR in $@; do
        [ -d $DIR -o -h $DIR -a $DIR ] && INFO "Found directory $DIR" || ERRO "NOT found directory $DIR"
    done
}

check_file() {
    for FILE in $@; do
        [ -f $FILE -o -h $FILE -a $FILE ] && INFO "Found file $FILE" || ERRO "NOT found file $FILE"
    done
}

dbfunc() {
    if [ $DEBUG -ne 0 ]; then
        echo -e $@
    else
        $@
    fi
}

# Here I use the new brand Hisat2, which could be computational efficient
# XXX: there could be a contradict when loading multiple tools
dbfunc source /apps/modules/modules.bashrc

# Project directories
proj_dir=~/Documents/projects/ASECausalSNPPrioritization
inpt_dir=$proj_dir/inputs
otpt_dir=$proj_dir/outputs
buf_dir=$proj_dir/buffers
temp_dir=$proj_dir/temps
check_dir $proj_dir $inpt_dir $otpt_dir $buf_dir $temp_dir

# The executable file
tlkt_dir=~/tools
hisat2_exe=$tlkt_dir/hisat2-2.1.0/hisat2
check_file $hisat2_exe

# Index
indx_dir=$buf_dir/hisat2_ss_exn_idx/human_g1k_v37

# Input
fstq_dir=$inpt_dir/BIOS_FASTQs

fstq_id=AC47H5ACXX-3-18
fstq_1=$fstq_dir/${fstq_id}_R1.fq.gz
fstq_2=$fstq_dir/${fstq_id}_R2.fq.gz
check_file $fstq_1 $fstq_2

# Output
bam_dir=$buf_dir/BIOS_alignments
temp_bam_dir=$temp_dir
sort_bam=$temp_bam_dir/${fstq_id}.sort.bam

hst2_thrd=10  # 10/15
smtl_thrd=5  # 5/15
dbfunc module load SAMtools
#
# Alignment using Hisat2 and transform the SAM output into BAM
# TODO: Should I mask loci (MAF>0.01) to remooth the reference bias?
#
dbfunc $hisat2_exe \
    -x $indx_dir \
    -1 $fstq_1 \
    -2 $fstq_2 \
    --rg-id $fstq_id \
    --rg 'SM:'$fstq_id \
    --rg 'PL:ILLUMINA' \
    --threads $hst2_thrd \
    | samtools sort \
    -@ $smtl_thrd \
    -O BAM \
    -o $sort_bam
check_file $sort_bam

sort_bai=${sort_bam/.bam/.bai}
dbfunc samtools index \
    -@ $smtl_thrd \
    $sort_bam \
    $sort_bai
check_file $sort_bai

#
# Some local calibration using GATK. Will it help?
#
dbfunc module load GATK/3.7-Java-1.8.0_74
mem=13  # 13/15
nt=15   # 15/15
nct=15  # 15/15

## Indel realignment
rfrn_genome=$inpt_dir/GRCh37_reference/human_g1k_v37.fasta
refg_dict=${rfrn_genome/.fasta/.dict}
knwn_indels=$inpt_dir/GATK_references/Mills_and_1000G_gold_standard.indels.b37.vcf
check_file $rfrn_genome $knwn_indels $refg_dict

# | --filter_reads_with_N_cigar / -filterRNC | Filter out reads with CIGAR containing the N operator, instead of failing with an error |
# | --threads / -nt                          | Number of file threads |
raln_itvl=${sort_bam/.bam/.raln.intervals}
dbfunc java -Xmx${mem}g -jar $EBROOTGATK/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -R $rfrn_genome \
    -I $sort_bam \
    -nt $nt \
    --known $knwn_indels \
    --filter_reads_with_N_cigar \
    -o $raln_itvl
check_file $raln_itvl

sort_raln_bam=${sort_bam/.bam/.raln.bam}
sort_raln_bai=${sort_bam/.bam/.raln.bai}
dbfunc java -Xmx${mem}g -jar $EBROOTGATK/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -R $rfrn_genome \
    -I $sort_bam \
    -o $sort_raln_bam \
    --knownAlleles $knwn_indels \
    --targetIntervals $raln_itvl \
    --filter_reads_with_N_cigar

check_file $sort_raln_bam
rm -fr $sort_bam $sort_bai $raln_itvl

## Base quality score recalibration
knwn_snps=$inpt_dir/GATK_references/dbsnp_138.b37.vcf
check_file $knwn_snps

bqsr_table=${sort_raln_bam/.bam/.bqsr.table}
dbfunc java -Xmx${mem}g -jar $EBROOTGATK/GenomeAnalysisTK.jar \
    -T BaseRecalibrator \
    -R $rfrn_genome \
    -I $sort_raln_bam \
    -o $bqsr_table \
    -nct $nct \
    -knownSites $knwn_snps
check_file $bqsr_table

sort_raln_bqsr_bam=${sort_raln_bam/.bam/.bqsr.bam}
sort_raln_bqsr_bai=${sort_raln_bam/.bam/.bqsr.bai}
dbfunc java -Xmx${mem}g -jar $EBROOTGATK/GenomeAnalysisTK.jar \
    -T PrintReads \
    -R $rfrn_genome \
    -I $sort_raln_bam \
    -o $sort_raln_bqsr_bam \
    -nct $nct \
    --BQSR $bqsr_table

check_file $sort_raln_bqsr_bam $sort_raln_bqsr_bai
rm -fr $sort_raln_bai $sort_raln_bam $bqsr_table

#
# Compile reads from the same read group into one BAM.
# But not applicable in current case
#

#
# Mark duplicates
#
dbfunc module load picard

sort_raln_bqsr_mkdp_bam=${sort_raln_bqsr_bam/.bam/.mkdp.bam}
sort_raln_bqsr_mkdp_txt=${sort_raln_bqsr_mkdp_bam/.bam/.dup_metrics.txt}

dbfunc java -Xmx${mem}g -jar $EBROOTPICARD/picard.jar \
    MarkDuplicates \
    I=$sort_raln_bqsr_bam \
    O=$sort_raln_bqsr_mkdp_bam \
    M=$sort_raln_bqsr_mkdp_txt

check_file $sort_raln_bqsr_mkdp_bam
rm -fr $sort_raln_bqsr_bai $sort_raln_bqsr_bam $sort_raln_bqsr_mkdp_txt

#
# Create an index file to 
#
dbfunc module load SAMtools

sort_raln_bqsr_mkdp_bai=${sort_raln_bqsr_mkdp_bam/.bam/.bai}
dbfunc samtools index \
    -@ $smtl_thrd \
    $sort_raln_bqsr_mkdp_bam \
    ${sort_raln_bqsr_mkdp_bai}

check_file $sort_raln_bqsr_mkdp_bai

mv $sort_raln_bqsr_mkdp_bam $sort_raln_bqsr_mkdp_bai $bam_dir

INFO "JOB was done"
