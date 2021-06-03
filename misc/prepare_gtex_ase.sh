#!/bin/bash

#
## A script to execute prepare_gtex_ase.py which generates ASE count table for given tissue.
#

# By default, this script submit jobs to the cluster workload manager SLURM, but it's also
# possible to run the script on local machine.

set -Eeu -o pipefail

# Workingspace
projdir=~/Documents/projects/wp_ase_dlp
iptdir=$projdir/inputs
optdir=$projdir/outputs
scpdir=$projdir/scripts

# Inputs
tissue=WHLBLD
ase_table_path=$iptdir/GTEx/ase
lookup_table=$iptdir/GTEx/misc/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.snp.only_id.tsv.gz

# Output
mkdir -p $optdir/aseQuan_v2/GTEx/$tissue

# local machine or cluster?
use_cluster=true
if [[ $use_cluster == true ]]; then
    # Resource
    reqmem=12G
    reqcpus=2
    reqtime=09:59:00
    jobname=prepare_gtex_ase
    logfile=$optdir/aseQuan_v2/GTEx/$tissue/%j-%u-prepare_gtex_ase.log

    execmd="sbatch \
        --time $reqtime \
        --mem-per-cpu $reqmem \
        --cpus-per-task $reqcpus \
        --output $logfile \
        --job-name $jobname"
else
    execmd='bash'
    reqcpus=2
fi

# Worker
$execmd <<EOF
#!/bin/bash
set -u -o pipefail

module load Python/3.7.4-GCCcore-8.3.0
source $scpdir/.env/bin/activate

# Determine the number of CPUs
threads=\${SLURM_CPUS_PER_TASK:-$reqcpus}

n_done_jobs=-\$threads
for ase_table in $ase_table_path/*.tsv.gz; do
    input_file=\$(basename \$ase_table)
    output_file=\${input_file/.tsv.gz/.$tissue.csv}

    # Skip current input file if the correponding output was found
    if [[ -f $optdir/aseQuan_v2/GTEx/$tissue/\$output_file ]]; then
        echo "[W]: Output file found, skip input file \$input_file"
        continue
    fi

    # If no given tissue was found, skip the ASE count table.
    has_tissue=\$(zgrep -wcm 1 $tissue \$ase_table)
    if [[ \$has_tissue -gt 0 ]]; then

        echo "[I]: Process input file \$input_file"
        python $scpdir/misc/prepare_gtex_ase.py \
            -l $lookup_table \
            -a \$ase_table \
            -t $tissue \
            -o $optdir/aseQuan_v2/GTEx/$tissue/\$output_file &

        n_done_jobs=\$(( n_done_jobs + 1 ))
    else
        echo "[W]: No record of $tissue was found in \$input_file."
    fi

    # Using the resource as efficient as possible.
    n_running_jobs=\$(jobs | grep -c Running)
    while [[ \$n_running_jobs -ge \$threads ]]; do
        echo -ne "[I]: \$n_running_jobs job(s) running. \$n_done_jobs job(s) done.\r"
        sleep 5s
        n_running_jobs=\$(jobs | grep -c Running)
    done
done

# Guard last background jobs.
if [[ \$n_running_jobs -gt 0 ]]; then wait; fi

echo "[I]: Totally, \$n_done_jobs job(s) were done."
EOF
