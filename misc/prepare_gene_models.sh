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

source ~/Documents/projects/wp_ase_dlp/scripts/.env/bin/activate

# Fetch Ensembl canonical and protein coding genes and transcripts.
python -c '#!/usr/bin/env python3
from argparse import ArgumentParser
import HTSeq

# NOTE: The GFF should be sorted.

parser = ArgumentParser()
parser.add_argument("-g", "--gff-file", required=True,
                    help="Input GFF file. Required")
parser.add_argument("-o", "--outfile", default="canonical_transcripts.gtf",
                    help="The output file name. Default: %(default)s")
opts = parser.parse_args()

gff = HTSeq.GFF_Reader(opts.gff_file, end_included=True)
cur_gene_rec = []
for idx, line in enumerate(gff):
    if line.type not in ["gene", "transcript", "exon"]: continue
    if line.attr["gene_type"] != "protein_coding": continue

    if line.type == "gene":
        cur_gene_rec.append(line)
    elif (hasattr(line, "attr") and "tag" in line.attr
        and "Ensembl_canonical" in line.attr["tag"]):
        cur_gene_rec.append(line)

with open(opts.outfile, "w") as outhandle:
    for line in cur_gene_rec:
        line = line.get_gff_line(with_equal_sign=True)
        line = line.replace("; ", ";").replace("\"", "").replace("chr", "")

        if line[0] not in ["X", "Y", "M"]:
            outhandle.write(line)
' -g $rf_dir/Gencode/gff3/gencode.v38.annotation.gff3.gz \
  -o $wk_dir/gencode.annotation.ensembl_canonical.protein_coding.GRCh38.gff3

# Liftover the obtained gene set coordinates from GRCh38 to GRCh37.
CrossMap.py gff $rf_dir/UCSCGenomeBrowser/chain/hg38ToHg19.over.nochr.chain \
  $wk_dir/gencode.annotation.ensembl_canonical.protein_coding.GRCh38.gff3 \
  $wk_dir/gencode.annotation.ensembl_canonical.protein_coding.GRCh37.gff3

# We also need to remove unmapped records
python -c '
import sys
bad_rec = []
with open(sys.argv[1]) as ipf:
    for line in ipf:
        rec_lst = line.split("\t")
        if len(rec_lst) != 9: continue
        for attr in rec_lst[-1].split(";"):
            if attr.startswith("gene_id"):
                bad_rec.append(attr.replace("gene_id=", ""))

with open(sys.argv[2]) as ipf, open(sys.argv[3], "w") as opf:
    for line in ipf:
        rec_lst = line.split("\t")
        if len(rec_lst) != 9: continue
        attr_dct = {k:v for k, v in [l.split("=") for l in rec_lst[-1].split(";")]}
        if attr_dct["gene_id"] in bad_rec: continue

        opf.write(line)
' $wk_dir/gencode.annotation.ensembl_canonical.protein_coding.GRCh37.gff3.unmap \
  $wk_dir/gencode.annotation.ensembl_canonical.protein_coding.GRCh37.gff3 \
  $wk_dir/gencode.annotation.ensembl_canonical.protein_coding.GRCh37.tmp.gff3

mv $wk_dir/gencode.annotation.ensembl_canonical.protein_coding.GRCh37.tmp.gff3 \
  $wk_dir/gencode.annotation.ensembl_canonical.protein_coding.GRCh37.gff3

# Remove the temporary file.
rm -f $wk_dir/gencode.v38.annotation.ensembl_canonical.protein_coding.gff3.gz
rm -f $wk_dir/gencode.annotation.ensembl_canonical.protein_coding.GRCh38.gff3
