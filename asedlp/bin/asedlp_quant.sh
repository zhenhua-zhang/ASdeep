#!/bin/bash
#SBATCH --time=0:29:0 \
#SBATCH --cpus=1 \
#SBATCH --mem=5G \
#SBATCH --array=1-22 \
#SBATCH --job-name=asedlp_quant \
#SBATCH --output=%A_%a-%u_asedlp_quant.log \

# Author     : Zhenhua Zhang
# Email      : zhenhua.zhang217@gmail.com
# License    : MIT
# Create date: Mon 09 Mar 2020 09:22:50 AM CET
# Last update: Mon 30 Mar 2020 04:31:45 PM CEST

set -o errexit
set -o errtrace


if [[ $(hostname) =~ "genetics" ]]; then
    echo "On local machine genetics."
else
    source /apps/modules/modules.bashrc
    module purge
    if [[ $(hostname) =~ "pg-node" || $(hostname) =~ "pg-gpu" ]]; then
        ml PyTorch/1.3.1-fosscuda-2019b-Python-3.7.4
    elif [[ $(hostname) =~ "umcg-node" ]]; then
        module load Python/3.6.3-foss-2015b
    fi
    module list
fi

# TODO: finish the help function
echo_help() {
    cat <<EOF

Usage:
    sbatch [sbatch-specific-options] $(basename "$0") [this-script-options]

Help:
    -w, --work-dir  Required.
        The work directory including outputs of rnaseqpl.

    -i, --fastq-id  Required.
        The the ID of fastq files.

    -a, --gff-file  Required.
        The genome feature format (GFF) file for reference genome. The file will
        be processed in to a SQLite3 database.

    -s, --sample-id-file  Required.
        The file including sample IDs: fastq_id-in-bios to sample-id-in-gonl.

    -G, --gene-id-file-dir  Required.
        The directory which includes files which gene IDs for each chrommosome.
        The file should be one line per gene ID and the gene IDs should be
        Ensembl ID which should also be consistent with ones in GFF files.

    -v, --venv-path  Optional.
        A Python virtual environment.

    -h, --help
        Print this help context and exit.

More information please contact Zhenhua Zhang <zhenhua.zhang217@gmail.com>

EOF
    exit 0
}

long_opts="work-dir:,fastq-id:,gff-file:,sample-id-file:,gene-id-file-dir:,venv:,help"
opt=$(getopt -l $long_opts -- "w:i:a:G:s:v:h" "$@")
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
		-h|--help) echo_help ;;
		--) shift && break ;;
	esac
	shift
done

work_dir=${work_dir:?-w/--work-dir is required!!}
fastq_id=${fastq_id:?-i/--fastq-id is required!!}
gff_file=${gff_file:?-g/--gff-file is required!!}
sample_id_file=${sample_id_file:?-s/--sample-id-file is required!!}
gene_id_file_dir=${gene_id_file_dir:?-G/--gene-id-file-dir is required!!}

# Ensure the job were executed under the job array mode of slurm>
if [ "${SLURM_ARRAY_TASK_ID}xxx" == "xxx" ]; then
    echo "[E]: the script should be run under job array mode."
    exit 1
fi

chrom_id=${SLURM_ARRAY_TASK_ID:=1}
sample_id=$(grep -w "$fastq_id" "$sample_id_file" | cut -f2)

if [ "${venv_path}xxx" != "xxx" ]; then
    source "$venv_path"
fi

ase_opt_dir=$work_dir/optdir/$fastq_id/aseOptdir
mkdir -p "$ase_opt_dir"/{ase_report,train_set}

gene_id_file=$gene_id_file_dir/chr$chrom_id.txt
haplotypes=$work_dir/snph5db/$chrom_id/haplotype.h5
snp_tab=$work_dir/snph5db/$chrom_id/snps_tab.h5
snp_index=$work_dir/snph5db/$chrom_id/snps_index.h5
sequence_tab=$work_dir/fastah5db/chr$chrom_id.h5
ref_read_counts=$work_dir/optdir/$fastq_id/waspOptdir/perChrom/$chrom_id/${fastq_id}_$chrom_id.refAlleleCounts.h5
alt_read_counts=$work_dir/optdir/$fastq_id/waspOptdir/perChrom/$chrom_id/${fastq_id}_$chrom_id.altAlleleCounts.h5
ase_report=$ase_opt_dir/ase_report/${fastq_id}_${chrom_id}_ase_report.txt
train_set=$ase_opt_dir/train_set/${fastq_id}_${chrom_id}_matrix_and_ase.npz

./asedlp quant \
    --sample-id "$sample_id" \
    --gene-id-file "$gene_id_file" \
    --haplotypes "$haplotypes" \
    --snp-tab "$snp_tab" \
    --snp-index "$snp_index" \
    --sequence-tab "$sequence_tab" \
    --ref-read-counts "$ref_read_counts" \
    --alt-read-counts "$alt_read_counts" \
    --genome-annot "$gff_file" \
    --save-as-ase-report "$ase_report" \
    --save-as-train-set "$train_set"
