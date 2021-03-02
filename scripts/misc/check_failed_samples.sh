#!/bin/bash

#
## Commands used to check the results.
#

# First check duplication. I found there are some duplicated results when
# checking optDir results.
projdir=~/Documents/projects/wp_ase_dlp

# LLS_660Q has many duplicates, after removing the duplicated ones,
# successful_count is equal to Schedualed_count now.
workdir=$projdir/outputs/aseQuan/logDir/LLS_660Q
for step in FastpPreproc STARMapping WaspRMBFindIntersectedSnps WaspRMBFilterAndRmDup WaspRMBRemapping CollectOutputAndCleanup; do
    while read -r item; do
            dup_files=( $(ls $workdir/*$item-$step.log.gz) );
            last=$(( ${#dup_files} - 1 ))
            file1_id=$(basename ${dup_files[0]} | tr '_' '-' | cut -f1 -d'-')
            file2_id=$(basename ${dup_files[$(( ${#dup_files[@]} - 1 ))]} | tr '_' '-' | cut -f1 -d'-')

            if [[ $file1_id -lt $file2_id ]]; then
                extra_file_id=$file1_id
            else
                extra_file_id=$file2_id
            fi

            case $step in
                WaspRMBFindIntersectedSnps | WaspRMBFilterAndRmDup )
                    mv $workdir/$extra_file_id"_"{1..22}* $(dirname $workdir)/LLS_660Q-duplicated 
                    continue ;;
                FastpPreproc | STARMapping | WaspRMBRemapping | CollectOutputAndCleanup)
                    mv $workdir/$extra_file_id-* $(dirname $workdir)/LLS_660Q-duplicated ;;
            esac
    done < <(ls $workdir/*$step.log.gz | xargs -n 1 basename | cut -f2-4 -d'-' | sort | uniq -d)
done

while read -r item; do
    count=$(grep -c $item"_" $projdir/scripts/configs/LLS_660Q-metarun.sh)
    if [[ $count -ne 0 ]]; then
        echo $item;
    fi
done < <(ls *.log.gz | cut -f2-4 -d'-' | sort | uniq)


#
## Table to check results
#
cat <<EOF
# Cohort_id	Successful_count	Schedualed_count	Pass
CODAM	180	180	Yes
LL	0	407	No
LLS_660Q	372	372	Yes
LLS_OminExpr	236	236	Yes
NTR_Affy6	744	744	Yes
PAN	167	167	Yes
RS	698	698	Yes
gonl	273	273	Yes
Total	2670	3077	No
EOF

# vim: set ts=25:
