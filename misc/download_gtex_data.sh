#!/bin/bash
projdir=~/Documents/projects/wp_ase_dlp

# The credentials.json has a expiration date, please check it.
~/tools/bin/gen3-client configure \
    --profile SWERTZ \
    --cred $projdir/inputs/GTExmisc/credentials.json \
    --apiendpoint=https://gen3.theanvil.io

for file_type in ase vcf; do
    ~/tools/bin/gen3-client download-multiple \
        --protocol s3 \
        --profile SWERTZ \
        --download-path $projdir/inputs/GTEx/$file_type \
        --manifest $projdir/inputs/GTEx/misc/$file_type-manifest.json
done
