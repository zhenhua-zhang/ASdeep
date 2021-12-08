# Author: Zhenhua Zhang
# E-mail: zhenhua.zhang217@gmail.com
# Created: May 07, 2020
# Updated: Jul 22, 2021

"""Predicting allele-specific expression by integrating RNA-seq read counts
and haplotypes.

Note:
    The package is an allele-specific expression analysis tool using the output
    from ASEReadCounter.
"""

import pdb
import argparse

from .inferai import inferai
from .database import makedb
from .predict import predict
from .train import train
from .zutils import LogManager
from .zutils import parse_verbose

def get_args():
    """Get CLI arguments."""
    parser = argparse.ArgumentParser(
        description="A tool to interpret variants function by deep-learning")
    parser.add_argument(
        "-l", dest="logfile",
        help="The file into which logging write.")
    parser.add_argument(
        "-v", dest="verbose", action="count", default=2,
        help="Verbose level. A counting keyword argument.")
    parser.add_argument(
        "-d", "--debug", action="store_true",
        help="Enter debug mode using pdb.")

    subparser = parser.add_subparsers(prog="asdeep", dest="subcmd")

    _parser = subparser.add_parser(
        "inferai", help="Infer allelic imbalance using allelic read counts.")
    _parser.add_argument(
        "-v", "--genetic-variants", metavar="FILE", required=True,
        help="Personal variants. Required")
    _parser.add_argument(
        "-i", "--genome-intervals", metavar="FILE", required=True,
        help="Gene or isoform position. Required")
    _parser.add_argument(
        "-c", "--readcounts-table", metavar="FILE", required=True,
        help="Metadata table in CSV format. Required")
    _parser.add_argument(
        "-s", "--sample-id", metavar="STR", required=True,
        help="The sample ID in the database to predict.")
    _parser.add_argument(
        "-d", "--n-draw", metavar="INT", default=1000, type=int,
        help="Number of draws in the MCMC. Default: %(default)s")
    _parser.add_argument(
        "-t", "--n-tune", metavar="INT", default=1000, type=int,
        help="Number of tunes in the MCMC. Default: %(default)s")
    _parser.add_argument(
        "-n", "--n-chain", metavar="INT", default=2, type=int,
        help="Number of MCMC chains. Default: %(default)s")
    _parser.add_argument(
        "-@", "--n-cpu", metavar="INT", default=2, type=int,
        help="Number of CPUs to burn. Default: %(default)s")
    _parser.add_argument(
        "-p", "--hdi-prob", metavar="FLOAT", default=0.90, type=float,
        help="The HDI to be claculated. Default: %(default)s")
    _parser.add_argument(
        "-f", "--feature", metavar="STR", default="exon",
        help="For which feature the allelic estimation. Default: %(default)s ")
    _parser.add_argument(
        "-o", "--out-file", metavar="PATH", default="allelic-expression.csv",
        help="Output directory. Default: %(default)s")

    _parser = subparser.add_parser(
        "makedb", help="Make a HDF5 database for train, test, and prediction.")
    _parser.add_argument(
        "-g", "--reference-genome", metavar="FILE", required=True,
        help="Reference genome. Required")
    _parser.add_argument(
        "-v", "--genetic-variants", metavar="FILE", required=True,
        help="Personal variants. Required")
    _parser.add_argument(
        "-i", "--genome-intervals", metavar="FILE", required=True,
        help="Gene or isoform position. Required")
    _parser.add_argument(
        "-m", "--metadata-table", metavar="FILE", required=True,
        help="Metadata table in CSV format. Required")
    _parser.add_argument(
        "-n", "--n-base-pairs", metavar="INT", default=262144, type=int,
        help="Number of base paires. Default: %(default)s")
    _parser.add_argument(
        "-o", "--out-dir", metavar="PATH", default="./",
        help="Output directory. Default: %(default)s")

    _parser = subparser.add_parser(
        "train", help="Train the model on quantified ASE effects")
    _parser.add_argument(
        "-d", "--database", required=True, metavar="FILE",
        help="The database including training dataset. Required")
    _parser.add_argument(
        "-a", "--prebuilt-arch", default="alexnet",
        choices=["alexnet", "resnext", "resnet"], metavar="STR",
        help="The prebuilt architecture to be used. Default: %(default)s")
    _parser.add_argument(
        "-m", "--model-state-path", default="asdeep_model.pth", metavar="FILE",
        help="The path to which save the model state. Default: %(default)s")
    _parser.add_argument(
        "-r", "--learning-rate", default=1e-5, type=float, metavar="FLT",
        help="Learning rate. Default: %(default)s")
    _parser.add_argument(
        "-e", "--epoches", default=5, type=int, metavar="INT",
        help="The number of epoch will be run. Default: %(default)s")
    _parser.add_argument(
        "-b", "--batch-size", default=64, type=int, metavar="INT",
        help="Batch size. Default: %(default)s")
    _parser.add_argument(
        "-n", "--n-base-pairs", default=65536, type=int, metavar="INT",
        help="Number of base pairs to be included in each training sample. "
        "Default: %(default)s")
    _parser.add_argument(
        "-p", "--train-pp", default=.9, type=float, metavar="FLT",
        help="The propotion of dataset will be used as training set. "
        "Default: %(default)s")
    _parser.add_argument(
        "-f", "--homo-flank", default=25,
        help="The flank length of heterozygous length. Default: %(default)s")
    _parser.add_argument(
        "-L", "--log-dir", default="./asdeep_logs", metavar="PATH",
        help="The path into which write the logging files for tensorboard. "
        "Default: %(default)s")
    _parser.add_argument(
        "-N", "--log-n-epoch", default=5, type=int, metavar="INT",
        help="Number of epoch for each evaluation point. Default: %(default)s")
    _parser.add_argument(
        "-R", "--random-state", default=42, type=int, metavar="INT",
        help="The random seed for torch. Default: %(default)s")
    _parser.add_argument(
        "-@", "--n-cpus", default=4, type=int, metavar="INT",
        help="Number workers for loading dataset. Default: %(default)s")

    _parser = subparser.add_parser(
        "predict", help="Predict based on given net work state")
    _parser.add_argument(
        "-d", "--database", required=True, metavar="FILE",
        help="The pattern to read input files.")
    _parser.add_argument(
        "-s", "--sample-ids", nargs="*", metavar="STR",
        help="The sample ID in the database to predict.")
    _parser.add_argument(
        "-a", "--prebuilt-arch", default="alexnet", metavar="STR",
        choices=["alexnet", "resnext", "resnet"],
        help="The prebuilt architecture to be used. Default: %(default)s")
    _parser.add_argument(
        "-m", "--model-state-path", default="asdeep_model.pth", metavar="FILE",
        help="The model state to be loaded. Default: %(default)s")
    _parser.add_argument(
        "-f", "--homo-flank", default=25, metavar="PATH",
        help="The flank length of heterozygous length. Default: %(default)s")
    _parser.add_argument(
        "-o", "--out-dir", default="asdeep_pred", metavar="PATH",
        help="The path to save the prediction results. Default: %(default)s")
    _parser.add_argument(
        "-t", "--attributions", nargs="*", metavar="STR",
        choices=["GS", "IG"],
        help="Whether show the attributions.")
    _parser.add_argument(
        "-F", "--save-fmt", default="png",
        help="The figure format for attributes. Default: %(default)s")

    return parser


def main():
    """Main function."""
    parser = get_args()
    args = parser.parse_args()

    logfile = args.logfile
    level = parse_verbose(args.verbose)
    logman = LogManager("Main", level=level, logfile=logfile)

    if args.debug:
        pdb.set_trace()

    if args.subcmd:
        logman.info("Job starts.")
        if args.subcmd == "inferai":
            inferai(args)
        elif args.subcmd == "makedb":
            makedb(args)
        elif args.subcmd == "train":
            train(args)
        elif args.subcmd == "predict":
            predict(args)
        logman.info("Job done!\n")
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
