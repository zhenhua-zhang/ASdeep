import os
import sys
import argparse

import numpy
import scipy
import pysam
import pandas
import matplotlib

from pysam import VariantFile
from pysam import AlignmentFile

def read_bam(bam_file):
    with open(bam_file, 'rb') as bam_file_handle:
        sam_context = AlignmentFile(bam_file_handle, "rb")
    

def read_bcf(bcf_file):
    with open(bcf_file, 'rb') as bcf_file_handle:
        vcf_context = VariantFile(bcf_file_handle)