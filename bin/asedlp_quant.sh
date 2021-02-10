#!/bin/bash
# Author     : Zhenhua Zhang
# Email      : zhenhua.zhang217@gmail.com
# License    : MIT
# Create date: Mon 09 Mar 2020 09:22:50 AM CET


cat <<EOF > /dev/null
pjdir=/groups/umcg-bios/tmp01/users/umcg-zzhang/projects/wp_ase_dlp
bash $pjdir/scripts/bin/asedlp_quant.sh \
    -i AC1JL5ACXX-1-10 \
    -w $pjdir/outputs/aseQuan \
    -a $pjdir/inputs/Ensembl_references/Homo_sapiens.GRCh37.75.gtf \
    -G $pjdir/inputs/Ensembl_references/gene_id-per-chr \
    -s $pjdir/inputs/idmapping/BIOS-genotype-rnaseq-ids_with-LL_usable_20210105.txt \
    -v $pjdir/scripts/.env
EOF

# Text color and format {
if [[ -x /usr/bin/tput ]] && tput setaf 1 &> /dev/null; then
    tput sgr0
    fmtBold=$(tput bold)
    fmtReset=$(tput sgr0)
    fmtItalic=$(tput sitm)

    fmtRed=$(tput setaf 1)
    fmtBlue=$(tput setaf 4)
    fmtCyan=$(tput setaf 6)
    fmtGrey=$(tput setaf 8)
    fmtBlack=$(tput setaf 0)
    fmtGreen=$(tput setaf 2)
    fmtWhite=$(tput setaf 7)
    fmtYellow=$(tput setaf 3)
    fmtPurple=$(tput setaf 5)
else
    fmtBold="\e[1m"
    fmtReset="\e[0m"
    fmtItalic="\e[3m"

    fmtRed="\e[1;31m"
    fmtBlue="\e[1;34m"
    fmtCyan="\e[1;36m"
    fmtGrey="\e[1;37m"
    fmtBlack="\e[1;30m"
    fmtGreen="\e[1;32m"
    fmtWhite="\e[1;37m"
    fmtYellow="\e[1;33m"
    fmtPurple="\e[1;35m"
fi
# } // Text color and format

# Error information, exit -1
echoErro() {
    echo -e "$fmtBold$fmtRed[E]: $1$fmtReset" >&2 && exit -1
}

# Warning information
echoWarn() {
    echo -e "$fmtBold$fmtYellow[W]:$fmtReset $*" >&2
}

# General information
echoInfo() {
    echo -e "$fmtBold$fmtWhite[I]:$fmtReset $*"
}

# Echo help for this script
echoHelp() {
    cat <<EOF

Usage:
    bash [sbatch-specific-options] $(basename "$0") [this-script-options]

Help:
    -w, --work-dir  Required.
        The work directory including outputs of rnaseqpl.

    -i, --fastq-id  Required.
        The the ID of fastq files.

    -a, --gff-file  Required.
        The genome feature format (GFF) file for reference genome. The file will
        be processed in to a SQLite3 database.

    -s, --sample-id-file  Required.
        The file including sample IDs: fastq_id-in-bios to sample-id-in-bios.

    -G, --gene-id-file-dir  Required.
        The directory which includes files which gene IDs for each chrommosome.
        The file should be one line per gene ID and the gene IDs should be
        Ensembl ID which should also be consistent with ones in GFF files.

    -v, --venv-path  Optional.
        A Python virtual environment.

    -d, --debug
        Debug mode

    -h, --help
        Print this help context and exit.

More information please contact Zhenhua Zhang <zhenhua.zhang217@gmail.com>

EOF
    exit 0
}

set -Ee

long_opts="work-dir:,fastq-id:,gff-file:,sample-id-file:,gene-id-file-dir:,venv-path:,debug,help"
opt=$(getopt -l $long_opts -- "w:i:a:G:s:v:dh" "$@")
eval set -- ${opt}
while true; do
	case $1 in
        # Required
        -w|--work-dir) shift && work_dir=$1 ;;
        -i|--fastq-id) shift && fastq_id=$1 ;;
        -a|--gff-file) shift && gff_file=$1 ;;
        -s|--sample-id-file) shift && sample_id_file=$1 ;;
        -G|--gene-id-file-dir) shift && gene_id_file_dir=$1 ;;

        # Optional
		-v|--venv-path) shift && venv_path=$1 ;;
        -d|--dubug) debug=1 ;;
		-h|--help) echoHelp ;;
		--) shift && break ;;
	esac
	shift
done

# Check options
fastq_id=${fastq_id:?-i/--fastq-id is required!!}
work_dir=${work_dir:?-w/--work-dir is required!!}
if [[ ! -d $work_dir ]]; then
    echoErro "Failed to find -w/--work-dir: $work_dir"
fi

gene_id_file_dir=${gene_id_file_dir:?-G/--gene-id-file-dir is required!!}
if [[ ! -d $gene_id_file_dir ]]; then
    echoErro "Failed to find -G/--gene-id-file-dir: $gene_id_file_dir"
fi

gff_file=${gff_file:?-a/--gff-file is required!!}
if [[ ! -e $gff_file ]]; then
    echoErro "Failed to find -a/--gff-file: $gff_file"
fi

sample_id_file=${sample_id_file:?-s/--sample-id-file is required!!}
if [[ ! -e $sample_id_file ]]; then
    echoErro "Failed to find -s/--sample-id-file: $sample_id_file"
fi


# # Ensure the job were executed under the job array mode of slurm
# if [ "${SLURM_ARRAY_TASK_ID}xxx" == "xxx" ]; then
#     echo "[E]: the script should be run under job array mode."
#     exit 1
# fi

# Config
chrom_id=${SLURM_ARRAY_TASK_ID:=1}
sample_id=$(grep -w "$fastq_id" "$sample_id_file" | cut -f2)
cohort_id=$(grep -w "$fastq_id" "$sample_id_file" | cut -f4)

# Input
gene_id_file=$gene_id_file_dir/chr$chrom_id.txt
seq_tab=$work_dir/fastaH5db/chr$chrom_id.h5
hap_tab=$work_dir/snpH5db/$cohort_id/$chrom_id/haplotype.h5
snp_tab=$work_dir/snpH5db/$cohort_id/$chrom_id/snps_tab.h5
snp_idx=$work_dir/snpH5db/$cohort_id/$chrom_id/snps_index.h5
ref_read_counts=$work_dir/optDir/$fastq_id/waspOptDir/perChrom/$chrom_id/${fastq_id}_$chrom_id.refAlleleCounts.h5
alt_read_counts=$work_dir/optDir/$fastq_id/waspOptDir/perChrom/$chrom_id/${fastq_id}_$chrom_id.altAlleleCounts.h5

# Output
ase_opt_dir=$work_dir/optDir/$fastq_id/aseOptDir
ase_report=$ase_opt_dir/aseReport/${fastq_id}_${chrom_id}_ase_report.txt
train_set=$ase_opt_dir/trainSet/${fastq_id}_${chrom_id}_matrix_and_ase.fa.gz

if [[ ${debug:=0} -eq 1 ]]; then
    cmd="bash"
else
    cmd="sbatch --time=0:19:0 --mem=5G --array=1-22 --cpus-per-task=1"
    cmd=$cmd" --job-name=${fastq_id}_asedlp_quant"
    cmd=$cmd" --output=$work_dir/logDir/%A_%a_${fastq_id}_asedlp_quant.log"
fi

$cmd <<EOF
#!/bin/bash

host_name=\$(hostname)
if [[ \${host_name} =~ "genetics" ]]; then
    echo "On local machine genetics."
else
    source /apps/modules/modules.bashrc
    module purge
    if [[ \${host_name} =~ "pg-node" || \${host_name} =~ "pg-gpu" ]]; then
        ml PyTorch/1.3.1-fosscuda-2019b-Python-3.7.4
    elif [[ \${host_name} =~ "gs-vcompute" || \${host_name} =~ "gearshift" ]]; then
        module load Python/3.7.4-GCCcore-7.3.0-bare
    elif [[ \${host_name} =~ "umcg-node" || \${host_name} =~ "calculon" ]]; then
        module load Python/3.6.3-foss-2015b
    fi
    module list
fi

if [ "${venv_path}xxx" != "xxx" ]; then
    source "$venv_path"/bin/activate
fi

mkdir -p "$ase_opt_dir"/{aseReport,trainSet}
~/Documents/projects/wp_ase_dlp/scripts/asedlp/asedlp quant \
    --sample-id "$sample_id" \
    --gene-id-file "$gene_id_file" \
    --haplotypes "$hap_tab" \
    --snp-tab "$snp_tab" \
    --snp-index "$snp_idx" \
    --sequence-tab "$seq_tab" \
    --ref-read-counts "$ref_read_counts" \
    --alt-read-counts "$alt_read_counts" \
    --genome-annot "$gff_file" \
    --save-as-ase-report "$ase_report" \
    --save-as-train-set "$train_set"
EOF
