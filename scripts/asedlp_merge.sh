#!/bin/bash

# File name : asedlp_merge.sh
# Author    : Zhenhua Zhang
# E-mail    : zhenhua.zhang217@gmail.com
# Created   : 2020 May 14
# Version   : v0.2.0
# License   : MIT

set -o errexit errtrace

pjdir=${HOME}/Documents/projects/wp_ase_dlp
ipdir=${pjdir}/inputs
opdir=${pjdir}/outputs
spdir=${pjdir}/scripts


sbatch --time 0:20:0 <<EOF
#!/bin/bash
#SBATCH --mem=2G
#SBATCH --time=0:20:0
#SBATCH --output=asedlp_merge.log
#SBATCH --job-name=asedlp_merge
#SBATCH --cpus-per-task=1


set -o errexit errtrace

# Check host name
case \$(hostname) in
    umcg-node* | calculon | boxy)
        source /apps/modules/modules.bashrc
        module purge && module load Python/3.6.3-foss-2015b && module list
        ;;
    pg-interactive* | pg-node* | pg-gpu*)
        module purge && module load PyTorch/1.3.1-fosscuda-2019b-Python-3.7.4 && module list
        ;;
    genetics)
        echo "On local machine genetics."
        ;;
    *)
        echo "[E]: Unknown hostname"
        exit
esac


if [ -d ${spdir}/.env ]; then
    source ${spdir}/.env/bin/activate
fi


${spdir}/asedlp/asedlp merge \
    -r \$(find -L ${pjdir}/workdir/optdir/ -name "*_ase_report.txt") \
    -p ${opdir}/ase_genes/individual_vs_gene_matrix/ase_genes \
    -g ${ipdir}/Ensembl_references/Homo_sapiens.GRCh37.75.gtf \
    --max-na-per-gene 80

for x in {1..22}; do
    ${spdir}/asedlp/asedlp merge \
        -r \$(find -L ${pjdir}/workdir/optdir/ -name "*_\${x}_ase_report.txt") \
        -p ${opdir}/ase_genes/individual_vs_gene_matrix/ase_genes.p_val_heatmap.chr\${x} \
        -g ${ipdir}/Ensembl_references/Homo_sapiens.GRCh37.75.gtf \
        --max-na-per-gene 50
done

EOF
