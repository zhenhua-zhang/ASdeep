#!/bin/bash
#SBATCH --time=:10:0
#SBATCH --mem=5G
#SBATCH --cpus=10
#SBATCH --output=%j-%u-fastp_preproc.log
#SBATCH --job-name=fastp_preproc

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


dbfunc source /apps/modules/modules.bashrc

# Project directories
proj_dir=~/Documents/projects/ASECausalSNPPrioritization
ipt_dir=$proj_dir/inputs
opt_dir=$proj_dir/outputs
buf_dir=$proj_dir/buffers
spt_dir=$proj_dir/scripts
check_dir $proj_dir $inpt_dir $otpt_dir $buf_dir $temp_dir

# The executable file
tlkt_dir=/groups/umcg-gcc/tmp03/umcg-zzhang/tools
fastp_exe=$tlkt_dir/fastp/fastp
check_file $fastp_exe

# Input FASTQ files
fastq_id=AC47H5ACXX-3-18
fastq_raw_dir=$ipt_dir/BIOS_FASTQs
fastq_1=$fastq_raw_dir/${fastq_id}_R1.fq.gz
fastq_2=$fastq_raw_dir/${fastq_id}_R2.fq.gz
check_file $fastq_1 $fastq_2

# Outputs. Paired reads, unpaired reads, failed reads, reports
fastq_clean_dir=$buf_dir/BIOS_FASTQs_clean
fastq_out_1=$fastq_clean_dir/${fastq_id}_paired_R1.fq.gz
fastq_out_2=$fastq_clean_dir/${fastq_id}_paired_R2.fq.gz
# fastq_unpaired_1=$fastq_clean_dir/${fastq_id}_unpaired_R1.fq.gz
# fastq_unpaired_2=$fastq_clean_dir/${fastq_id}_unpaired_R2.fq.gz
# fastq_failed=$fastq_clean_dir/${fastq_id}_failed.fq.gz
fastq_report=$fastq_clean_dir/${fastq_id}_report.html
check_dir $fastq_clean_dir

# Other settings
fastp_threads=10

#
# Preprocessing reads using fastp
# --detect_adapter_for_pe is enabled to check the adapters, it could slow down the process by the manual.
# --trim_poly_g is enabled to trim poly-G, where the G means no signal in the Illumina two-color systems.
# --overrepresentation_sampling is used to check 1% of all reads to analysis the overrepresented sequence
#
$fastp_exe \
    --in1 $fastq_1 \
    --in2 $fastq_2 \
    --out1 $fastq_out_1 \
    --out2 $fastq_out_2 \
    --html $fastq_report \
    --detect_adapter_for_pe \
    --cut_front \
    --cut_tail \
    --correction \
    --trim_poly_g \
    --trim_poly_x \
    --overrepresentation_sampling 100 \
    --thread $fastp_threads \
    --dont_overwrite

    # --unpaired1 $fastq_unpaired_1 \
    # --unpaired2 $fastq_unpaired_2  \
    # --failed_out $fastq_failed \