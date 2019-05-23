#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Get ASE effects for each sample"""
import itertools
import os

from argparse import ArgumentParser

import HTSeq

from utilities import print_arguments


TRANS_MATRIX = {
    "AA":(1, 0, 0, 0), "AC":(1, 1, 0, 0), "AG":(1, 0, 1, 0), "AT":(1, 0, 0, 1),
    "aa":(1, 0, 0, 0), "ac":(1, 1, 0, 0), "ag":(1, 0, 1, 0), "at":(1, 0, 0, 1),
    "CA":(1, 1, 0, 0), "CC":(0, 1, 0, 0), "CG":(0, 1, 1, 0), "CT":(0, 1, 0, 1),
    "ca":(1, 1, 0, 0), "cc":(0, 1, 0, 0), "cg":(0, 1, 1, 0), "ct":(0, 1, 0, 1),
    "GA":(1, 0, 0, 1), "GC":(0, 1, 1, 0), "GG":(0, 0, 1, 0), "GT":(0, 0, 1, 1),
    "ga":(1, 0, 0, 1), "gc":(0, 1, 1, 0), "gg":(0, 0, 1, 0), "gt":(0, 0, 1, 1),
    "TA":(1, 0, 0, 1), "TC":(0, 1, 0, 1), "TG":(0, 0, 1, 1), "TT":(0, 0, 0, 1),
    "ta":(1, 0, 0, 1), "tc":(0, 1, 0, 1), "tg":(0, 0, 1, 1), "tt":(0, 0, 0, 1),
    "NN":(0, 0, 0, 0), "nn":(0, 0, 0, 0)
}


def get_args():
    """A method to get arguments from the command line"""

    parser = ArgumentParser()
    _group = parser.add_argument_group("Global")
    _group.add_argument(
        "-V", dest="verbose_level", action="count", help="Verbose level"
    )

    subparser = parser.add_subparsers(dest="subcommand")
    _group = parser.add_argument_group("Input")
    _group.add_argument(
        "--input-file", dest="input_file", type=str, required=True,
        help="Input file. Required"
    )

    return parser

class ASEReadCounter(HTSeq):
    """ASE reads counter"""

    def __init__(self, sam_file, gff_file):
        """Init of ASEReadCounter"""
        super(ASEReadCounter, self).__init__(file_name)
        self.sam_hand = self.read_bam(bam_file)
        self.gff_hand = self.read_gff(gff_file)


class SequenceMatrix(FileReader):
    """Make sequence matrix"""

    def __init__(self, file_name):
        """Init of SequenceMatrix"""
        super(SequenceMatrix, self).__init__(file_name)


# Construction of parental genomes
class ParentalGenome():
    """Construct parental genomes from .vcf file and .fasta file"""

    def __init__(self, variants, sequences):
        """Init"""
        self.variants = variants
        self.sequences = sequences

    def parse_variants(self):
        """Parse variant file"""

    def parse_sequences(self):
        """Parse sequence file"""

    def make_parental_genome(self):
        """Construct parental genomes"""

