#!/bin/bash
# Author : Zhenhua Zhang
# E-mail : zhenhua.zhang217@gmail.com
# Created: Sep 10, 2021
# Updated: Sep 10, 2021

# Preapre a gene set in GFF format.

# The transcripts are:
#   1. Ensembl Canonical (tagged by Ensembl_canonical) 
#   2. Protein coding (type is protein_coding)
#   3. Under the GRCh38 coordinates but then, liftovered to GRCh37 coordinates.
#   4. all from autosomes.

# TODO: The gene set should include exon features.
set -Eeu -o pipefail

# Reference and working directories.
rf_dir=~/Documents/projects/wp_reference
wk_dir=~/Documents/projects/wp_ase_dlp/inputs/Gencode

# Fetch Ensembl canonical and protein coding genes and transcripts.
in_gff=$rf_dir/Gencode/gtf/gencode.v38.annotation.gtf.gz 
sub_gtf=$wk_dir/gencode.annotation.ensembl_canonical.protein_coding.GRCh38.gtf
v37_gtf=$wk_dir/gencode.annotation.ensembl_canonical.protein_coding.GRCh37.gtf
chain_file=$rf_dir/UCSCGenomeBrowser/chain/hg38ToHg19.over.nochr.chain

awk -F$'\t' -f- <(zcat $in_gff) <<'EOF' | sed 's/^chr//' > $sub_gtf
/##/ {print; next}
{
  if ($1 ~ /chr[MXY]/) { next }

  if ($9 ~ /\<protein_coding\>/)
  {
    if ($3 == "gene") { print }
    else if ($3 == "transcript" || $3 == "exon") {
      if ($9 ~ /\<Ensembl_canonical\>/) {print}
    }
  }
}
EOF

# Liftover the obtained gene set coordinates from GRCh38 to GRCh37.
CrossMap.py gff $chain_file $sub_gtf $v37_gtf

# We also need to remove unmapped records
# NOTE: A small trick to sort gene/transcript/exon, check key_array line.
awk -F$'\t' -f- $v37_gtf*.unmap $v37_gtf <<'EOF' \
  | sort -t$'\t' -k2,2 -k5,5n -k1,1 -k6,6n \
  | cut -f2- -d$'\t' \
  >| $(dirname $v37_gtf)/tmp.gtf
BEGIN {key_array["gene"]=1; key_array["transcript"]=2; key_array["exon"]=3}
FNR == NR {
  if ($3 == "gene")
  {
    split($9, attr_dict, ";")
    split(attr_dict[1], gene_id, " ")
    ensg_id = gensub(/"(.+)"/, "\\1", "g", gene_id[2])

    unmap[ensg_id] = 1
  }

  next
}
{
  split($9, attr_dict, ";")
  split(attr_dict[1], gene_id, " ")
  ensg_id = gensub(/"(.+)"/, "\\1", "g", gene_id[2])

  if (unmap[ensg_id] != 1) {print key_array[$3]"\t"$0}
}
EOF

mv -f $(dirname $v37_gtf)/tmp.gtf $v37_gtf

# Compress and index
module purge
module load HTSlib/1.10.2-GCCcore-7.3.0

bgzip $v37_gtf
tabix -f -p gff $v37_gtf.gz

# Remove the temporary file.
rm -f $sub_gtf $v37_gtf*.unmap
