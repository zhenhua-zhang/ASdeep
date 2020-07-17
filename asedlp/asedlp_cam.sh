#!/bin/bash
#
# File name : asedlp_cam.sh
# Author    : Zhenhua Zhang
# E-mail    : zhenhua.zhang217@gmail.com
# Created   : Mon 22 Jun 2020 10:46:41 AM CEST
# Version   : 800
# License   : MIT
#

set -eu -o pipefail

pjdir=${HOME}/Documents/projects/wp_ase_dlp
ipdir=${pjdir}/inputs
opdir=${pjdir}/outputs
spdir=${pjdir}/scripts

model_state_dir=${opdir}/models/version_0.1.0

# file_pat=${pjdir}/'workdir/optdir/**/aseOptdir/train_set/*_1_matrix_and_ase.npz'
# gene_id=ENSG00000163221
#
# file_pat=${pjdir}/'workdir/optdir/**/aseOptdir/train_set/*_2_matrix_and_ase.npz'
# gene_id=ENSG00000173272
#
# file_pat=${pjdir}/'workdir/optdir/**/aseOptdir/train_set/*_15_matrix_and_ase.npz'
# gene_id=ENSG00000116710
#
# file_pat=${pjdir}/'workdir/optdir/**/aseOptdir/train_set/*_16_matrix_and_ase.npz'
# gene_id=ENSG00000206177

file_pat=${pjdir}/'workdir/optdir/**/aseOptdir/train_set/*_17_matrix_and_ase.npz'
# gene_id=ENSG00000108405
# gene_id=ENSG00000269871
# gene_id=ENSG00000170315
# gene_id=ENSG00000185862
gene_id=ENSG00000161570


sbatch --time 0:19:0 <<EOF
#!/bin/bash
#SBATCH --mem=5G
#SBATCH --ntasks=1
#SBATCH --output=${sbatch_log_path}
#SBATCH --job-name=${gene_id}.asedlp_train


set -eu -o pipefail

# Check host name
case \$(hostname) in
    umcg-node* | calculon | boxy)
        source /apps/modules/modules.bashrc
        module purge
        module load Python/3.6.3-foss-2015b
        module list
        ;;
    pg-interactive* | pg-node* | pg-gpu*)
        module purge
        module load PyTorch/1.3.1-fosscuda-2019b-Python-3.7.4
        module load TensorFlow/2.1.0-fosscuda-2019b-Python-3.7.4
        module list
        ;;
    genetics)
        echo "On local machine genetics."
        ;;
    *)
        echo "[E]: Unknown hostname"
        exit
esac


${spdir}/asedlp/asedlp cam \
    --gene-id ${gene_id} \
    --model-state ${model_state_path} \
    --file-path ${file_pat}

EOF
