#!/usr/bin/env python
'''
This is a utility script to match genotypes and RNA-seq results for BIOS
project.
'''
import sys

rnaseqfile = sys.argv[1]
with open(rnaseqfile, 'r') as fh:
    rnaseqs = list({x.rsplit("/", 1)[-1].split("_")[0] for x in fh.readlines()})

idmappingfile = sys.argv[2]
with open(idmappingfile, 'r') as fh:
    ids = [x.strip().split("\t") for x in fh.readlines()]
    ids_dict = {dna_id: rna_id.split("_", 1)[0] for dna_id, rna_id in ids}
    assert len(ids) == len(ids_dict)


genotypes = {}
genotypefile = sys.argv[3]
with open(genotypefile, 'r') as fh:
    for item in fh.readlines():
        biobk_id, dna_id = item.strip().split("\t")
        if biobk_id == "CODAM":
            key = "-".join([biobk_id, dna_id.split("_", 1)[-1]])
        elif biobk_id == "gonl":
            key = dna_id
        elif biobk_id == "LL":
            key = "-".join([biobk_id, dna_id.split("_", 1)[-1]])
        elif biobk_id in ["LLS_660Q", "LLS_OminExpr"]:
            key = "-".join([biobk_id.split("_")[0], dna_id.split("_", 1)[-1]])
        elif biobk_id == "NTR_Affy6":
            key = "-".join([biobk_id.split("_")[0], dna_id.split("_", 1)[-1]])
        elif biobk_id == "PAN":
            key = "-".join([biobk_id, dna_id.split("_", 1)[-1]])
        elif biobk_id == "RS":
            key = "-".join([biobk_id, dna_id.split("_", 1)[-1]])

        genotypes[key] = [dna_id, biobk_id]

match_ids = [[key, dna_id, ids_dict[key], bionk_id]
             for key, (dna_id, bionk_id) in genotypes.items()
             if ids_dict.get(key, "NULL") in rnaseqs]

print("#uniq_id\tdna_id\trna_id\tbiobank_id")
for item in match_ids:
    print("\t".join(item))
