#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author      : Zhenhua Zhang
# Email       : zhenhua.zhang217@gmail.com
# License     : MIT
# Create date : 2020-03-07
# Last update : 2020-03-07

import argparse

from ased.ASEFactory import ASEFactory
from ased.DLPFactory import CNNFactory


def get_args():
    """Get CLI arguments.
    """
    parser = argparse.ArgumentParser()
    sub_parsers = parser.add_subparsers(prog="asedlp", dest="subcmd", description="A factory to produce ASE effects.")

    quant_parser = sub_parsers.add_parser("quant", help="Quantification of ASE effects")
    quant_parser.add_argument("-i", "--sample-id", action="store", default="gonl-101b", dest="sample_id", help="The sample ID.")
    quant_parser.add_argument("-g", "--gene-ids", action="store", nargs="+", default="ENSG00000008130", dest="gene_ids", help="The gene IDs.")
    quant_parser.add_argument("-G", "--gene-id-file", action="store", default=None, dest="gene_id_file", help="The file from which read the gene ids.")
    quant_parser.add_argument("--haplotypes", action="store", default="./haps.h5", dest="haplotypes", help="Path to HDF5 file to read phased haplotypes from.")
    quant_parser.add_argument("--snp-index", action="store", default="./snp_idex.h5", dest="snp_index", help="Path to HDF5 file to read SNP index from.")
    quant_parser.add_argument("--snp-tab", action="store", default="./snp_tab.h5", dest="snp_tab", help="Path to HDF5 file to read SNP information from.")
    quant_parser.add_argument("--sequence-tab", action="store", default="./seq_tab.h5", dest="seq_tab", help="Path to HDF5 file to read reference sequence from.")
    quant_parser.add_argument("--ref-read-counts", action="store", default="./ref_count_tab.h5", dest="ref_tab", help="Path to HDF5 file to read reference reads counts from.")
    quant_parser.add_argument("--alt-read-counts", action="store", default="./alt_count_tab.h5", dest="alt_tab", help="Path to HDF5 file to read alternative reads counts from.")
    quant_parser.add_argument("--genome-annot", action="store", default="./genome.gff.gz", dest="genome_annot", help="Path to GFF / GTF file to read gene structure information from.")

    train_parser = sub_parsers.add_parser("train", help="Train the model on quantfified ASE effects.")
    train_parser.add_argument("-i", "--input-file", action="store", default="./trainset.npz", dest="train_set", help="The file from which read the training dataset.")

    return parser.parse_args()

def quant(args):
    factory = ASEFactory(args)
    factory.init() \
            .gen_gnm_region() \
            .gen_seq_mtrx() \
            .gen_ase_effect() \
            .save_ase_report() \
            .save_for_training() \
            .shutdown()

def train(args):
    pass

def main():
    args = get_args()
    if args.subcmd == "quant":
        quant(args)
    elif args.subcmd == 'train':
        train(args)
    else:
        raise Exception("Unknown subcommand")

if __name__ == '__main__':
    main()
