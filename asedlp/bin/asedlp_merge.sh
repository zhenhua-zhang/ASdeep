#!/bin/bash
#SBATCH --mem=2G
#SBATCH --time=0:20:0
#SBATCH --output=asedlp_merge.log
#SBATCH --job-name=asedlp_merge
#SBATCH --cpus-per-task=1

set -o errexit
set -o errtrace

# Check host name
host_name=$(hostname)
if [[ $host_name =~ "genetics" ]]; then
    echo "On local machine genetics."
else
    module purge
    case $host_name in
        umcg-node* | calculon | boxy)
            source /apps/modules/modules.bashrc
            module load Python/3.6.3-foss-2015b
            ;;
        pg-interactive* | pg-node* | pg-gpu*)
            module load Python/3.7.2-GCCcore-8.2.0
            ;;
        *)
            echo "[E]: Unknown hostname"
            exit
    esac
    module list
fi

if [ -d ../../.env ]; then
    source ../../.env/bin/activate
fi


pjdir=~/Documents/projects/wp_ase_dlp
opdir=${pjdir}/outputs/ase_genes/individual_vs_gene_matrix

./asedlp merge \
    -r $(find -L ${pjdir}/workdir/optdir/ -name "*_ase_report.txt") \
    -p ${opdir}/ase_genes.p_val_heatmap \
    -g ../../../inputs/Ensembl_references/Homo_sapiens.GRCh37.75.gtf \
    --max-na-per-gene 50

for x in {1..22}; do
    ./asedlp merge \
        -r $(find -L ${pjdir}/workdir/optdir/ -name "*_${x}_ase_report.txt") \
        -p ${opdir}/ase_genes.p_val_heatmap.chr${x} \
        -g ../../../inputs/Ensembl_references/Homo_sapiens.GRCh37.75.gtf \
        --max-na-per-gene 50
done
