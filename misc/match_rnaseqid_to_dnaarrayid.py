#!/usr/bin/env python
'''
getids
'''
import sys

rec_order = [
    "ids", "bios_id", "uuid", "biobank_id.x", "person_id", "pheno_id",
    "biobank_gwas_id", "dna_id", "rna_id.x", "rna_note", "gonl_id",
    "old_gonl_id", "cg_id", "in_rp3", "biobank_id.y", "rna_id.y", "run_id",
    "type", "freeze", "qc", "raw", "past_filter"
]

def mk_conv_dict(conv_rec, header_idx=0, sep='\t'):
    header = conv_rec.pop(header_idx).strip().split(sep)
    conv_dict = {}
    for rec in conv_rec:
        reclist = rec.strip().split(sep)
        val = dict(zip(header, reclist))
        if rec.startswith('NTR'):
            key = 'NTR-' + reclist[0].split('-')[-1]
        else:
            key = reclist[0]

        if key not in conv_dict:
            conv_dict[key] = val
        else:
            print(key)
            raise KeyError("Duplicated keys")

    return conv_dict


def mk_search_key(dnaid_rec, sep="\t"):
    searchkeys = {}
    for rec in dnaid_rec:
        rec = rec.strip()
        gntpid = rec.split(sep)[-1]
        if rec.startswith("CODAM"):
            key = rec.replace('\t', '-')
        elif rec.startswith("NTR"):
            key = "NTR-" + rec.split("F")[-1]
        elif rec.startswith("PAN"):
            key = "PAN-" + rec.split("_")[-1]
        elif rec.startswith("RS"):
            key = "RS-" + rec.split("_")[-1]
        else:
            raise ValueError("Unknown record: {}".format(rec))

        searchkeys[key] = gntpid

    return searchkeys


def mt_ids(cdict, skeys):
    for key in skeys:
        if key in cdict:
            _convrec = cdict[key]
            print(skeys[key] + "\t" + "\t".join(_convrec[x] for x in rec_order))
        else:
            print("Missing key-val pair: {}, {}".format(key, skeys[key]))


if len(sys.argv) != 3:
    sys.exit('Wrong number of arguments: require two, given {}'.format(len(sys.argv)))

bios_id_conv_file = sys.argv[1]     # BIOS ID convert table
bios_dnaarray_id_file = sys.argv[2] # BIOS genotype id

with open(bios_id_conv_file) as convfh, open(bios_dnaarray_id_file, 'r') as dnafh:
    convrec = convfh.readlines()
    dna_array_rec = dnafh.readlines()

    convdict = mk_conv_dict(convrec, sep=",")
    searchkey = mk_search_key(dna_array_rec)
    mt_ids(convdict, searchkey)
