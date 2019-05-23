#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Make matrix"""
# gonl-100a AC47H5ACXX-3-18 Example
# 1. get some example data from real dataset
# 1.1 example.fa
# reference genome: GRCh37
# (/groups/umcg-bios/tmp03/users/umcg-zzhang/projects/ASEPrediction/benchmark/inputs/references)
# ```
# $> head -668 genome.fa > tmp.fa
# $> grep -vn N tmp.fa
# $> (head -1 tmp.fa; sed -n '169,$p' tmp.fa) > example.fa
# $> module load SAMtools/1.5-foss-2015b
# $> samtools faidx example.fa
# example.fa example.fa.fai
# ```
# Now we have a fragment reference genome from GRCh37. It starts with non-N at
# 10,021 (167 * 60), ends at 40,020 (667 * 60), with length 30,000. (1:1-40020)
# 1.2 example.bam
# ```
# module load SAMtools
# ```
# 1.3 example.vcf.gz
# Fetch variants from gonl.chr1.snps_index.r5.3.vcf.gz
# /groups/umcg-gonl/prm02/releases/variants/GoNL1/release6.1/06_IL_haplotype_panel
# ```
# $> module load BCFtools/1.6.foss-2015b
# $> bcftools view gonl.chr1.snps_indels.r5.3.vcf.gz 1:40021-70020 > ~/Documents/example.vcf
# $> bgzip ~/Documents/example.vcf
# $> cd ~/Documents
# $> bcftools index example.vcf.gz
# ```

# 1.4 example.gff.gz
# Fetch interval from Homo_sapiens.GRCh37.75.gff
# ```
# grep -P "^1\t" Homo_sapiens.GRCh37.75.gff | sort -k4,5n > example.gff
# module load tabix
# bgzip example.gff
# tabix -p gff example.gff.gz
# awk '{if ($5 <= 9250800 && $1 == 1) {print}}' ../gff/Homo_sapiens.GRCh37.75.gtf  > example.gff
# ```

# Sequence is 0-based


import os
import sys

import pysam

fas_file = "test/example.fa"
vcf_file = "test/example.vcf.gz"
gtf_file = "test/example.gtf.gz"
bam_file = "test/example.bam"

class UnknownStrandError(Exception):
    pass

class LengthError(Exception):
    pass

class NonDictError(Exception):
    pass

# A made matrix transforming alleles into matrix, which will make life easier
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

def parse_variant(variants):  # TODO: check data type of var.pos and var.chr
    record_hash = {
        (var.chrom, var.pos) : [
            var.chrom, var.pos, var.id, var.ref, var.alts, var.qual, var.info
        ] for var in variants
    }  # FIXME: could be duplicated positions
    return record_hash

def matrix_yielder(seq, var_hash, chrom, shift):
    # FIXME: `shift` should be determined by strand???
    if isinstance(var_hash, dict):
        var_keys = var_hash.keys()
    else:
        raise NonDictError  # TODO: a concrete sub-class of TypeError

    for _idx, base in enumerate(seq):
        pos = _idx + shift

        if (chrom, pos) in var_keys:
            record = var_hash[(chrom, pos)]  # TODO: func to handle multi alts
            base = record[3] + record[4][0]
        else:
            base += base

        if len(base) == 2:
            yield TRANS_MATRIX[base]
        else:  # TODO: Pipe stderr into a log file
            print(
                "Non-single base, will use reference for both alleles: "
                "{}".format(base), file=sys.stderr
            )

def make_matrix(inter_hand, seq_hand, var_hand, with_orf=False, contig="1",
                up_shift=1000, dw_shift=1000, merged=True):
    """Make input matrix"""
    var_header = var_hand.header
    # var_sample = var_hand.samples

    for field in inter_hand.fetch("1"):
        # Important fields
        attrbutes = field.attributes
        gene_id = field.gene_id
        contig = field.contig
        strand = field.strand
        start = field.start
        end = field.end

        up_start = (int(start) - 1) - up_shift
        up_end = int(start) - 1

        dw_start = int(end) + 1
        dw_end = (int(end) + 1) + dw_shift

        up_seq = seq_hand.fetch(reference=contig, start=up_start, end=up_end)
        up_var = var_hand.fetch(
            contig=contig, start=up_start, stop=up_end, reopen=True
        )
        up_var_hash = parse_variant(up_var)
        up_matrix = matrix_yielder(up_seq, up_var_hash, contig, up_start)

        dw_seq = seq_hand.fetch(reference=contig, start=dw_start, end=dw_end)
        dw_var = var_hand.fetch(
            contig=contig, start=dw_start, stop=dw_end, reopen=True
        )
        dw_var_hash = parse_variant(dw_var)
        dw_matrix = matrix_yielder(dw_seq, dw_var_hash, contig, dw_start)

        # TODO: should yield a merged metrix of upstream and downstream
        if strand == '-':
            yield (dw_matrix, up_matrix)

        yield (up_matrix, dw_matrix)


sequence_hand = pysam.FastaFile(fas_file)
variant_hand = pysam.VariantFile(vcf_file, duplicate_filehandle=True)
alignment_hand = pysam.AlignmentFile(bam_file)
interval_hand = pysam.TabixFile(gtf_file, parser=pysam.asGTF())

matrix_pool = make_matrix(interval_hand, sequence_hand, variant_hand)

for up_dw_matrix in matrix_pool:
    up_matrix, dw_matrix = up_dw_matrix
    for matrix in up_matrix:
        print(list(matrix))

    for matrix in dw_matrix:
        print(list(matrix))

    break

sequence_hand.close()
variant_hand.close()
alignment_hand.close()
interval_hand.close()
