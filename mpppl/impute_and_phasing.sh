#!/bin/bash
#
# File name : impute_and_phasing.sh
# Author    : Zhenhua Zhang
# E-mail    : zhenhua.zhang217@gmail.com
# Created   : Tue 30 Jun 2020 04:49:32 PM CEST
# Version   : V 0.1.0
# License   : MIT
#

sbatch --time 1:59:0 \
    --mem 16G \
    --cpus-per-task 15 \
    --output %j-%u-shapeit_chr1.log \
    --job-name shapeit <<EOF
#!/bin/bash

set -eu -o pipefail

source /apps/modules/modules.bashrc
pjdir=\${HOME}/Documents/projects/wp_ase_dlp
ipdir=\${pjdir}/inputs
opdir=\${pjdir}/outputs
bfdir=\${pjdir}/buffers

module load shapeit
shapeit -B chr1 \
    -M \${ipdir}/1kg_p3/gmap/chr1.b37.gmap.gz \
    -T 15 \
    -O \${bfdir}/chr1.phased \
    --seed 31415

shapeit -convert \
    --input-haps \${bfdir}/chr1.phased \
    --output-ref \${bfdir}/chr1.phased.hap \${bfdir}/chr1.phased.leg \${bfdir}/chr1.phased.sam
EOF

sbatch --time 0:20:0 \
    --mem 16G \
    --cpus-per-task 8 \
    --output %j-%u-impute2_chr1.log \
    --job-name shapeit <<EOF
#!/bin/bash

set -eu -o pipefail

source /apps/modules/modules.bashrc
pjdir=\${HOME}/Documents/projects/wp_ase_dlp
ipdir=\${pjdir}/inputs
opdir=\${pjdir}/outputs

module load IMPUTE2
impute2 \
    -use_prephased_g \
    -m \${ipdir}/1kg_p3/1000GP_Phase3/genetic_map_chr1.txt \
    -h \${ipdir}/1kg_p3/1000GP_Phase3/1000GP_Phase3_chr1.hap.gz \
    -l \${ipdir}/1kg_p3/1000GP_Phase3/1000GP_Phase3_chr1.legend.gz \
    -known_haps_g \${pjdir}/buffers/chr1.phased.haps \
    -Ne 20000 \
    -o chr1.phased.impute2 \
    -int 1 5e6 \
    -phase

    # -strand_g Strand alignments it not usually necessary when you just want to
    # phase a dataset, but it is important when that dataset will be combined
    # with a reference panel in a down stram analysis
EOF
