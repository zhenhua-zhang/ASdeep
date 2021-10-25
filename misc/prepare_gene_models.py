#!/usr/bin/env python3
from argparse import ArgumentParser
import HTSeq

# NOTE:
#   1. The GFF should be sorted

def main():
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
            outhandle.write(line.get_gff_line(with_equal_sign=True)
                            .replace("; ", ";").replace("\"", ""))


if __name__ == "__main__":
    main()
