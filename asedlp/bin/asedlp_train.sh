#!/bin/bash
#SBATCH --mem=5G
#SBATCH --time=0:59:0
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --partition=gpu
#SBATCH --output=%j.%u.batch_asedlp_train.log
#SBATCH --job-name=batch_asedlp_train
#
# File name : asedlp_train.sh
# Author    : zhzhang
# E-mail    : zhzhang2015@sina.com
# Created   : Wed 01 Apr 2020 04:28:18 PM CEST
# Version   : v0.0.1
# License   : MIT

# Run `asedlp train` which is a python script using job array by slurm.
# Example:
# pjdir=~/Documents/projects/ASEDLO

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

file_pat=${file_pat:='../../../workdir/optdir/**/aseOptdir/train_set/*_17_matrix_and_ase.npz'}
gene_id=${gene_id:=ENSG00000108405}
model_state_path=${model_state_path:=${gene_id}.ptz}

./asedlp train \
    --file-pat "${file_pat}" \
    --gene-id "${gene_id}" \
    --cv-times "${cv_times:=8}" \
    --n-epoch "${n_epoch:=15}" \
    --model-state-path "${model_state_path}"
