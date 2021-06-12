#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: May 07, 2020
# Updated: Jun 07, 2021

"""Predicting allele-specific expression by integrating RNA-seq read counts and haplotypes.

Note:
    The package is a allele-specific expression analysis tool using the output from
    ASEReadCounter.

TDOO:
"""

import sys
import logging
import argparse


def get_args():
    """Get CLI arguments."""
    parser = argparse.ArgumentParser(description="A tool for the analysis of allele-specific expression.")
    parser.add_argument("-l", "--log-file", dest="log_file", default=None, metavar="FILE",
                        help="The file into which logging write.")
    parser.add_argument("-v", "--verbose-level", action="count", dest="verbose_level", default=2,
                        help="Verbose level. A counting keyword argument.")

    sub_parsers = parser.add_subparsers(prog="asedeep", dest="subcmd")

    quant_parser = sub_parsers.add_parser("quantify", help="Quantify ASE effects and generate training set")
    quant_parser.add_argument("-i", "--sample-id", required=True, dest="sample_id", metavar="STRING",
                              help="The sample ID in variants file.")
    quant_parser.add_argument("-v", "--variants-from", required=True, dest="variants_file", metavar="FILE",
                              help="The file from which read variants information.")
    quant_parser.add_argument("-s", "--genome-seq-from", required=True, dest="genome_seq_file", metavar="FILE",
                              help="The file from which read the genome sequence.")
    quant_parser.add_argument("-f", "--gene-feature-from", required=True, dest="gene_feature_file", metavar="FILE",
                              help="The file (GFF/GTF) from which read the gene feature.")
    quant_parser.add_argument("-r", "--ase-readcount-from", required=True, dest="ase_readcount_file", metavar="FILE",
                              help="The file from which read the allele-specific read counts.")
    quant_parser.add_argument("-g", "--gene-ids", default=[], nargs="+", dest="gene_ids", metavar="STRING",
                              help="The gene IDs. Default: %(default)s")
    quant_parser.add_argument("-G", "--gene-id-from", default=None, dest="gene_id_file", metavar="FILE",
                              help="The file from which read the gene ids. Default: %(default)s")
    quant_parser.add_argument("-m", "--min-allelic-counts", default=3, dest="min_allelic_counts", metavar="INTEGER",
                              help="The minimal read counts from each chromosome. Default: %(default)s")
    quant_parser.add_argument("-M", "--min-total-counts", default=6, dest="min_total_counts", metavar="INTEGER",
                              help="The minimal read counts from two chromosomes. Default: %(default)s")
    quant_parser.add_argument("--seq-only", action="store_true", dest="seq_only",
                              help="Disable sequence generation.")
    quant_parser.add_argument("--ase-only", action="store_true", dest="ase_only",
                              help="Disable ASE generation.")
    quant_parser.add_argument("-T", "--train-set-to", default="ase-train-set.fa", dest="as_train_set", metavar="FILE",
                              help="The file to which write the training dataset into. Default: %(default)s")
    quant_parser.add_argument("-R", "--ase-report-to", default="ase-report.txt", dest="as_ase_report", metavar="FILE",
                              help="The file to which write the ASE effect into. Default: %(default)s")
    quant_parser.add_argument("-F", "--shift-factor", default=5e5, type=int, dest="shift_factor", metavar="INTEGER",
                              help="Number of base pairs shift from the start codon, the strand will be considered. Default: %(default)s")

    report_parser = sub_parsers.add_parser("report", help="Report the quantified ASE effects")
    report_parser.add_argument("-p", "--file-path", required=True, dest="file_path_pool", nargs="+", metavar="FILE",
                               help="The pattern to read input files.")
    report_parser.add_argument("-f", "--gene-feature-db", required=True, dest="gene_feature_db", metavar="FILE",
                               help="The file (GFF/GTF) from which read the gene feature.")
    report_parser.add_argument("-e", "--exclude-from", default=None, dest="exclude_from", metavar="FILE",
                               help="Exclude gene ID from the file.")
    report_parser.add_argument("-o", "--save-path", default="asedeep_report", dest="save_path", metavar="PATH",
                               help="The path into which save the prediction results. Default: %(default)s")

    train_parser = sub_parsers.add_parser("train", help="Train the model on quantified ASE effects")
    train_parser.add_argument("-p", "--file-path", required=True, dest="file_path_pool", metavar="FILE", nargs="+",
                              help="The pattern to read input files.")
    train_parser.add_argument("-g", "--gene-id", required=True, nargs="+", dest="gene_id_pool", metavar="STRING",
                              help="The Ensembl id of a gene for which train the model. Default: %(default)s")
    train_parser.add_argument("-s", "--model-state-path", default="./asedeep_model.pth", dest="model_state_path", metavar="FILE",
                              help="The path to which save the model state dict. Default: %(default)s")
    train_parser.add_argument("-L", "--logging-path", default="./asedeep_logs", dest="logging_path", metavar="DIRECTORY",
                              help="The path into which write the logging files for tensorboard. Default: %(default)s")
    train_parser.add_argument("-r", "--learning-rate", default=1e-5, type=float, dest="learning_rate", metavar="FLOAT",
                              help="Learning rate. Default: %(default)s")
    train_parser.add_argument("-e", "--n-epoch", default=5, type=int, dest="n_epoch", metavar="INTEGER",
                              help="The number of epoch will be run. Default: %(default)s")
    train_parser.add_argument("-R", "--random-state", default=42, type=int, dest="random_state", metavar="INTEGER",
                              help="The random seed for torch. Default: %(default)s")
    train_parser.add_argument("-b", "--batch-size", default=64, type=int, dest="batch_size", metavar="INTEGER",
                              help="Batch size. Default: %(default)s")
    train_parser.add_argument("-a", "--prebuilt-arch", default="alexnet", dest="prebuilt_arch", choices=["alexnet", "resnext", "resnet"], metavar="STRING",
                              help="The prebuilt architecture to be used. Default: %(default)s")
    train_parser.add_argument("--pp-train", default=.9, type=float, dest="pp_train", metavar="FLOAT",
                              help="The propotion of dataset will be used as training set. Default: %(default)s")
    train_parser.add_argument("--log-per-n-epoch", default=5, type=int, dest="log_per_n_epoch", metavar="INTEGER",
                              help="Number of epoch for each evaluation point. Default: %(default)s")
    train_parser.add_argument("--loss-curve-path", default="./", type=str, dest="loss_curve_path", metavar="PATH",
                              help="The path to save training loss curve. Default: %(default)s")
    train_parser.add_argument("-n", "--n-cpus", default=12, type=int, dest="n_cpus", metavar="INTEGER",
                              help="Number workers for loading dataset. Default: %(default)s")

    pred_parser = sub_parsers.add_parser("predict", help="Predict based on given net work state")
    pred_parser.add_argument("-p", "--file-path", required=True, dest="file_path_pool", nargs="*", metavar="FILE",
                             help="The pattern to read input files.")
    pred_parser.add_argument("-g", "--gene-id", default="all", dest="gene_id", metavar="STRING",
                             help="The Ensembl id of a gene for which train the model. Default: %(default)s")
    pred_parser.add_argument("-s", "--model-state-path", default="./asedeep_model.pth", dest="model_state_path", metavar="FILE",
                             help="The model state to be loaded. Default: %(default)s")
    pred_parser.add_argument("-a", "--prebuilt-arch", default="alexnet", dest="prebuilt_arch", choices=["alexnet", "resnext", "resnet"], metavar="STRING",
                             help="The prebuilt architecture to be used. Default: %(default)s")
    pred_parser.add_argument("-o", "--save-path", default="asedeep_pred", dest="save_path", metavar="PATH",
                             help="The path into which save the prediction results. Default: %(default)s")
    pred_parser.add_argument("-A", "--show-attributions", action="store_true", dest="show_attr",
                             help="Whether show the attributions.")

    return parser


def pickup_model(prebuilt_arch):
    """Pick up a modified models."""
    try:
        import torch.nn as nn
        import torchvision.models as models
    except ImportError as ime:
        logging.error("Import error: {}".format(ime))
        sys.exit(1)

    if prebuilt_arch == "resnext":
        # By default, using the pre-built ResNext.
        net = models.resnext50_32x4d(pretrained=False)
        net.conv1 = nn.Conv2d(1, 64, 11, 2, 3, bias=False)
        net.fc = nn.Linear(2048, 3, bias=True)
    elif prebuilt_arch == "resnet":
        net = models.resnet18(pretrained=False)
        net.conv1 = nn.Conv2d(1, 64, 11, 2, 3, bias=False)
        net.fc = nn.Linear(512, 3, bias=True)
    else:
        if prebuilt_arch != "alexnet":
            logging.warning("Unsupported prebuilt architecture, using default AlexNet!")
        # Now I moved to AlexNet
        net = models.alexnet(pretrained=False)
        net.features[0] = nn.Conv2d(1, 64, 11, 4, 2)
        net.classifier[6] = nn.Linear(4096, 3, bias=True)

    return net


def quantify(args):
    """Qunatify the allele-specific expression effect for genes."""
    try:
        from asedeep.quantify import Quantifier
    except ImportError as ime:
        logging.error("Import error: {}".format(ime))
        sys.exit(1)

    seq_only, ase_only = args.seq_only, args.ase_only
    if seq_only and ase_only:
        raise ValueError("--seq-only and --ase-only are mutually exclusive!")
    elif seq_only:
        Quantifier(args).gen_gnm_itvl().gen_seq().save_train_set()
    elif ase_only:
        Quantifier(args).gen_gnm_itvl().gen_ase().save_ase_report()
    else:
        Quantifier(args).gen_gnm_itvl().gen_seq().gen_ase().save_ase_report().save_train_set()


def report(args):
    """Generate a report for the quantified ASE effects."""
    try:
        from asedeep.report import ASEReport
    except ImportError as ime:
        logging.error("Import error: {}".format(ime))
        sys.exit(1)

    save_path = args.save_path
    exclude_from = args.exclude_from
    file_path_pool = args.file_path_pool
    gene_feature_db = args.gene_feature_db

    black_list = []
    if exclude_from is not None:
        with open(exclude_from) as ifhandle:
            black_list = [x.strip() for x in ifhandle]

    (ASEReport(file_path_pool, save_path, gene_feature_db, black_list=black_list)
     .save_pval_matrices()
     .save_figures(figheight=4, figwidth=24))


def train(args):
    """Train a CNN model."""
    try:
        import numpy
        import torch
        from asedeep.train import Trainer
    except ImportError as ime:
        logging.error("Import error: {}".format(ime))
        sys.exit(1)

    file_path_pool = args.file_path_pool
    gene_id_pool = args.gene_id_pool
    n_cpus = args.n_cpus
    n_epoch = args.n_epoch
    pp_train = args.pp_train
    batch_size = args.batch_size
    random_state = args.random_state
    logging_path = args.logging_path
    learning_rate = args.learning_rate
    prebuilt_arch = args.prebuilt_arch
    log_per_n_epoch = args.log_per_n_epoch
    loss_curve_path = args.loss_curve_path
    model_state_path = args.model_state_path

    numpy.random.seed(random_state)
    torch.manual_seed(random_state)
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True

    net = pickup_model(prebuilt_arch)

    logging.info("Batch size:    {}".format(batch_size))
    logging.info("Random state:  {}".format(random_state))
    logging.info("Epoch number:  {}".format(n_epoch))
    logging.info("Architecture:  {}".format(prebuilt_arch))
    logging.info("Training prop: {}".format(pp_train))
    logging.info("Learning rate: {}".format(learning_rate))

    (Trainer(net, gene_id_pool, file_path_pool, logging_path, log_per_n_epoch, n_cpus)
     .init()
     .load_dataset()
     .train_test_split(pp_train)
     .train(eps=n_epoch, learning_rate=learning_rate, batch_size=batch_size)
     .loss_curve(loss_curve_path)
     .save_model(model_state_path))


def predict(args):
    """Predict."""
    try:
        from asedeep.predict import Predictor
    except ImportError as ime:
        logging.error("Import error: {}".format(ime))
        sys.exit(1)

    file_path_pool = args.file_path_pool
    gene_id = args.gene_id
    save_path = args.save_path
    show_attr = args.show_attr
    prebuilt_arch = args.prebuilt_arch
    model_state_path = args.model_state_path

    net = pickup_model(prebuilt_arch)

    (Predictor(gene_id=gene_id, file_path=file_path_pool, net=net)
     .init(model_state=model_state_path)
     .load_dataset()
     .predict(save_path, show_attr=show_attr))


def main():
    """Main function."""
    parser = get_args()
    args = parser.parse_args()

    log_file = args.log_file
    verbose_level = args.verbose_level * 10
    logging.basicConfig(filename=log_file,
                        format="{levelname: ^8}| {asctime} | {message}",
                        style="{",
                        datefmt="%Y%m%d,%H:%M:%S",
                        level=verbose_level)

    if args.subcmd == "quantify":
        quantify(args)
    elif args.subcmd == "report":
        report(args)
    elif args.subcmd == "train":
        train(args)
    elif args.subcmd == "predict":
        predict(args)
    else:
        parser.print_help()
