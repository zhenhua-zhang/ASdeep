#!/bin/bash
#SBATCH --time=23:59:00
#SBATCH --mem=5G
#SBATCH --cpus-per-task=2
#SBATCH --output=download_igsr_data.log
#SBATCH --job-name=download_igsr_data

# As there was a user storage limitation (500G) at '/data/umcg-asedeep', I switch to a scheme of
# 'download then process'.

set -Eeu
wget_dw() {
    local help_str list_file url_col mdt_col url md5 opts

    help_str="wget_dw -f LIST_FILE [-u COLNUM] [-m COLNUM] [-l LOG_FILE] [-O OUTPUT_DIR] [-h]"
    opts=$(getopt -l "list-file:,url-col:,md5-col:,log-file:,output-dir,help" -- "f:u:m:l:O:h" "$@")
    eval set -- $opts
    while true; do
        case $1 in
            -f|--list-file)
                shift && list_file=$1 ;;
            -l|--log-file)
                shift && log_file=$1 ;;
            -u|--url-col)
                shift && url_col=$1 ;;
            -m|--md5-col)
                shift && md5_col=$1 ;;
            -O|--output-dir)
                shift && out_dir=$1 ;;
            -h|--help)
                echo -e "$help_str" && exit -1 ;;
            --)
                shift && break ;;
        esac
        shift
    done

    list_file=${list_file:?-l/--list-file is required!}
    log_file=${log_file:=wget_download.log}
    out_dir=${out_dir:=./}
    url_col=${url_col:=0}
    md5_col=${md5_col:=1}

    while read -a item; do
        file_url=${item[$url_col]}
        file_md5=${item[$md5_col]}

        if [[ $file_url =~ ^'#' ]]; then
            continue
        fi

        if [[ $file_url =~ ftp || $file_url =~ http ]]; then
            fname=$(basename $file_url)
            wget $file_url -O $out_dir/$fname -a $out_dir/$log_file
            real_md5=$(md5sum $out_dir/$fname | cut -f1 -d' ')
            [[ "$file_md5" == "$real_md5" ]] \
                && echo -e "Pass    $real_md5  $file_md5  $fname" \
                || echo -e "Failed  $real_md5  $file_md5  $fname"
        fi
    done < "$list_file"
}

pjdir=~/Documents/projects/wp_ase_dlp
dw_genotypes=false
dw_rnaseq=true

#
## Download Geuvadis genotypes
#
if [[ $dw_genotypes == true ]]; then
    wget_dw -f $pjdir/inputs/geuvadis/samples/igsr_Phase3_integrated-variant-call-sets.tsv \
        -O $pjdir/inputs/geuvadis/genotypes
fi

#
## Download fastq files: ERR188
#
if [[ $dw_rnaseq == true ]]; then
    wget_dw -f $pjdir/inputs/geuvadis/samples/igsr_Geuvadis-mRNA_ERR188.tsv \
        -O $pjdir/inputs/geuvadis/fastqs
fi

# vim: set ai nowrap ts=4 tw=200:
