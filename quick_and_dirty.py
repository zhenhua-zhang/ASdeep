import argparse
import gzip
import io
import os
import sys

import matplotlib
import numpy
import pandas
import pysam
import scipy
from pysam import AlignmentFile, TabixFile, VariantFile

_bgzf_magic = b"\x1f\x8b\x08\x04"

base_dict = dict(
    A='A', T='T', C='C', G='G', X='AT', Y='AC', Z='AG', O='TC', P='TG', Q='CG'
)


class FileIOError(Exception):
    def __init__(self, msg=None):
        if msg is None:
            msg = "An error occured for file i/o"
        super(FileIOError, self).__init__(msg)

class FileNotReadableError(Exception):
    def __init__(self, msg=None):
        if msg is None:
            msg = "File not readable!!"
        super(FileNotReadableError, self).__init__(msg)

class FileTypeError(Exception):
    def __init__(self, msg=None):
        if msg is None:
            msg = "Unsupported file type"
        super(FileTypeError, self).__init__(msg)

def isSAM(file_name):
    return file_name.endswith('.sam')

def isBAM(file_name):
    return file_name.endswith('.bam')

def isCRAM(file_name):
    return file_name.endswith('.cram')

def isVCF(file_name):
    return file_name.endswith('.vcf')

def isVCFGZ(file_name):
    return file_name.endswith('.vcf.gz')

def isBCF(file_name):
    return file_name.endswith('.bcf')

def isFASTA(file_name):
    return file_name.endswith('.fa') or file_name.endswith('.fasta')

def isFASTAGZ(file_name):
    return file_name.endswith('.fa.gz') or file_name.endswith('.fasta.gz')

def isGFF(file_name):
    return file_name.endswith('.gff')

def isGFFGZ(file_name):
    return file_name.endswith('.gff.gz')

def isGTF(file_name):
    return file_name.endswith('.gtf')

def isGTFGZ(file_name):
    return file_name.endswith('.gtf.gz')

def isBED(file_name):
    return file_name.endswith('.bed')

def isBEDGZ(file_name):
    return file_name.endswith('.bed.gz')



def read_interval(fn):
    if isGFF(fn):
        pass
    elif isGFFGZ(fn):
        pass
    elif isGTF(fn):
        pass
    elif isGTFGZ(fn):
        pass
    elif isBED(fn):
        pass
    elif isBEDGZ(fn):
        pass

    with TabixFile(fn) as file_handle:
        return file_handle


def check_index(fn):
    if isVCFGZ(fn):
        idex_ext = ".tbi"
    elif  isFASTAGZ(fn):
        idex_ext = '.fai'
    elif isBAM(fn):
        idex_ext = ".bai"
    elif isCRAM(fn):
        idex_ext = '.crai'
    else:
        raise FileTypeError("Unknow file type: {}".format(fn))

    fnidx = fn +idex_ext

    if not os.access(fnidx, os.F_OK):
        return False

    if not os.access(fnidx, os.R_OK):
        raise FileNotReadableError(
            "Found file but NOT readable: {}".format(fnidx)
        )

    return True

def create_preset(fn):
    if isGFFGZ(fn):
        return "gff"
    elif isBED(fn):
        return "bed"
    elif isSAM(fn):
        return "sam"
    elif isVCFGZ(fn):
        return "vcf"
    else:
        raise FileIOError("Uknown type of files: {}".format(fn))


def make_tabidx(fn):
    global _bgzf_magic
    if fn.endswith('.gz'):
        with open(fn, 'rb') as _bgzf:
            if _bgzf_magic == _bgzf.read(4):
                print(
                    "NOT a BGZF file, will recompress it and try again",
                    file=sys.stderr
                )
            else:
                # TODO: function to recompress the file, but keep the original
                pass
    else:  # Compress the file but keep the umpressed one
        try:
            pysam.tabix_compress(fn, fn + '.gz', force=True)
            fn += '.gz'
        except Exception as e:
            raise e

    try:
        prst = create_preset(fn)
        pysam.tabix_index(fn, force=True, preset=prst)
    except Exception as e:
        raise e

def read_variant(fn):
    if isBCF(fn):
        mode = 'rb'
    elif isVCF(fn) or isVCFGZ(fn):
        mode = 'r'

    if check_index(fn) == False:
        make_tabidx(fn)

    try:
        variants = VariantFile(fn, mode)
    except Exception as exc:
        raise exc

    return variants


def read_alignment(fn):
    if isBAM(fn):
        mode = 'rb'
    elif isCRAM(fn):
        mode = 'rc'

    try:
        alignments = AlignmentFile(fn, mode)
    except Exception as exc:
        raise exc

    return alignments


def count_reads(alignment_file, interval, min_depth, min_quality, avr_quality):
    pass

def make_sequence_matrix(sequence):
    pass

def likelihood_func():
    pass

def cal_likelihood(ase_counts):
    likelihood = 100
    return likelihood

def main(align_fn, var_fn, seq_fn, inte_fn, target_sample, configs):
    variants = read_variant(variant_fn)
    variants_header = variants.header
    variants_contigs = variants_header.contigs
    variants_samples = variants_header.samples
