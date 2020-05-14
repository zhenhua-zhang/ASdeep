#!/bin/bash
#SBATCH --mem=5G
#SBATCH --time=1:29:0
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --partition=gpu
#SBATCH --output=/home/p281736/Documents/projects/wp_ase_dlp/outputs/models/version_0.1.0/%j.%u.batch_asedlp_train.log
#SBATCH --job-name=batch_asedlp_train

# File name : asedlp_train.sh
# Author    : Zhenhua Zhang
# E-mail    : zhenhua.zhang217@gmail.com
# Created   : Wed 01 Apr 2020 04:28:18 PM CEST
# Version   : v0.0.1
# License   : MIT


set -o errexit
set -o errtrace

# TODO: finish the help function
echo_help() {
    cat <<EOF

Usage:

Help:

    -h, --help
        Print this help context and exit.

More information please contact Zhenhua Zhang <zhenhua.zhang217@gmail.com>

EOF
    exit 0
}

version="0.1.0"
while [[ "$1" =~ ^- && ! "$1" == "--" ]]; do
    case $1 in
        -V | --version )
            echo $version && exit ;;
        -p | --file-pat )
            shift && file_pat=$1 ;;
        -g | --gene-id )
            shift && gene_id=$1 ;;
        -s | --model-state-path)
            shift && model_state_path=$1 ;;
        -c | --cv-times)
            shift && cv_times=$1 ;;
        -n | --n-epoch)
            shift && n_epoch=$1 ;;
        -h | --help)
            echo_help ;;
        --)
            shift && break ;;
    esac
    shift
done

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
            module load PyTorch/1.3.1-fosscuda-2019b-Python-3.7.4
            ;;
        *)
            echo "[E]: Unknown hostname"
            exit
    esac
    module list
fi

model_state_dir=/home/p281736/Documents/projects/wp_ase_dlp/outputs/models/version_0.1.0

file_pat='../../../workdir/optdir/**/aseOptdir/train_set/*_17_matrix_and_ase.npz'
gene_id=ENSG00000108405

file_pat='../../../workdir/optdir/**/aseOptdir/train_set/*_2_matrix_and_ase.npz'
gene_id=ENSG00000173272
model_state_path=${model_state_path:=${model_state_dir}/${gene_id}.model_state.ptz}
loss_curve_path=${loss_curve_path:=${model_state_dir}/${gene_id}}

./asedlp train \
    --file-pat "${file_pat}" \
    --gene-id "${gene_id}" \
    --n-epoch "${n_epoch:=400}" \
    --model-state-path "${model_state_path}" \
    --loss-curve-path "${loss_curve_path}"
