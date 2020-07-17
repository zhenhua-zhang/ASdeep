#!/bin/bash
#
# File name : asedlp_train.sh
# Author    : Zhenhua Zhang
# E-mail    : zhenhua.zhang217@gmail.com
# Created   : 2020 Apr 01
# Version   : v0.2.0
# License   : MIT

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
gene_id=ENSG00000108405
gene_id=ENSG00000269871
gene_id=ENSG00000170315
gene_id=ENSG00000161570
# gene_id=ENSG00000185862

output_dir=${model_state_dir}/${gene_id}
mkdir -p "${output_dir}"

n_epoch=50
time_stamp=$(date +%Y%m%d_%H%M%S)
model_state_path=${output_dir}/${time_stamp}.model_state.ptz
logging_path=${output_dir}/${time_stamp}.tensorboard_log
sbatch_log_path=${output_dir}/${time_stamp}.%j.%u.asedlp_train.log

# loss_curve_path=${loss_curve_path:=${model_state_dir}/${gene_id}}


sbatch --time 0:29:0 <<EOF
#!/bin/bash
#SBATCH --mem=5G
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --partition=gpu
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


${spdir}/asedlp/asedlp train \
    --file-path ${file_pat} \
    --gene-id ${gene_id} \
    --n-epoch ${n_epoch} \
    --model-state-path ${model_state_path} \
    --logging-path ${logging_path} \
    --log-per-n-epoch 1

EOF

# --loss-curve-path ${loss_curve_path} \
