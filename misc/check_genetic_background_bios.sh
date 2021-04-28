#!/bin/bash
#
# File name : genetic_variants_pca.sh
# Author    : Zhenhua Zhang
# E-mail    : zhenhua.zhang217@gmail.com
# Created   : 2021-04-21
# Version   : V0.1.0
# License   : MIT
#

set -Eeu -o pipefail

if [[ -e /apps/modules/modules.bashrc ]]; then
    source /apps/modules/modules.bashrc
fi

# Cohorts: CODAM  NTR  PAN  RS
# Sample information
cohort="${1:?Error: You should give the cohort id!}"

# Sub dirs
#pjdir=/Documents/projects/wp_ase_dlp
pjdir=/groups/umcg-bios/tmp03/users/umcg-zzhang/projects/wp_ase_dlp
if [[ ! -d $pjdir/inputs/BIOS_genotypes/$cohort ]]; then
    echo "Failed to find $pjdir/inputs/BIOS_genotypes/$cohort! Exit"
    exit
fi
