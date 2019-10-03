#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Make matrix"""
import os
import sys

import numpy as np
import pysam as ps
import pandas as pd


def get_args():
    """Get CLI arguments for current script"""
    parser = ArgumentParser()
    _group = parser.add_argument_group("Input")
    _group.add_argument("-r", "--reference", type=str, dest="reference", action="store_value", required=True, help="The reference genome in FASTA format. Required")
    _group.add_argument("-v", "--variants", type=str, dest="variants", action="store_value", required=True, help="The vairants file in VCF format. Required")
    _group.add_argument("-a", "--annotations", type=str, dest="annotations", action="store_value", required=True, help="The annotations file in GFF / GTF format. Required")
    _group.add_argument("-s", "--alignments", type=str, dest="alignments", action="store_value", required=True, help="The sequence alignments file in SAM / BAM format. Required")
    # _group = parser.add_argument_group("Config")
    # _group = parser.add_argument_group("Output")
    return parser

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

NT2VEC = { # Encoded neucliotide
    "A": [1, 0, 0, 0], "C": [0, 1, 0, 0], "G": [0, 0, 1, 0], "T": [0, 0, 0, 1]
}

VEC2NT = {
    (1, 0, 0, 0): "A", (0, 1, 0, 0): "C", (0, 0, 1, 0): "G", (0, 0, 0, 1): "T"
}

DNTVEC = [0, 0, 0, 0]  # Base not in ACGT

NT2AMB = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "AC": "K", "AG": "Y", "AT": "W", "CG": "S", "CT": "R", "GT": "M",
    "CA": "k", "GA": "y", "TA": "w", "GC": "s", "TC": "r", "TG": "m",
}

AMB2NT = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "K": "AC", "Y": "AG", "W": "AT", "S": "CG", "R": "CT", "M": "GT",
    "k": "CA", "y": "GA", "w": "TA", "s": "GC", "r": "TC", "m": "TG",
}

def _make_variant_dict(variants):
    record_dict = {(var.chrom, var.pos): var for var in variants}  # Using all filter: PASS / Inaccessibl. FIXME: could be duplicated positionse
    return record_dict

def _parse_sample_genotype(samples, alleles):
    sample_dict = {key: _decode_vcf_haps(alleles, val.get("GT")) for key, val in samples.items()}
    return sample_dict

def _decode_vcf_haps(alleles, code):
    return (alleles[code[0]], alleles[code[1]])

def _encode_hap_into_vec(hap):
    return NT2VEC.get(hap, [0, 0, 0, 0])

def _decode_vec_into_hap(vec, dft="N"):
    return VEC2NT.get(tuple(vec), dft)

def matrix_factory(seq, var_hash, chrom, shift, target_samples=["gonl-100a"]):
    # FIXME: `shift` should be determined by strand???
    if not isinstance(var_hash, dict):
        raise NonDictError  # TODO: a concrete sub-class of TypeError

    if isinstance(target_samples, str):
        target_samples = [target_samples]

    target_sample_allele_vec = {}
    for _each_sample in target_samples:
        for (_chrom, _pos), _vcf_rec in var_hash.items():
            _ref_allele = seq[_pos]  # Reference allele base in sequence
            _alleles = _vcf_rec.alleles  # tuple of reference allele followed by alt alleles
            assert _ref_allele == _alleles[0]

            _phase = _vcf_rec.samples.get(_each_sample, None).get("GT", None)
            if _phase[0] == _phase[1]:
                continue

            seq = seq[:_pos] + NT2AMB[_alleles[_phase[0]], _alleles[_phase[1]]]] + seq[_pos+1:]

        _allele_vec = []
        for _base in seq:
            _allele_0, _allele_1 = AMB2NT[_base]
            _allele_vec.append(NT2VEC.get(_allele_0, DNTVEC) + NT2VEC.get(_allele_1, DNTVEC]))

        target_sample_allele_vec[_each_sample] = _allele_vec
        print(seq)
    
    return target_sample_allele_vec

def make_matrix(inter_hand, seq_hand, var_hand, aln_hand, with_orf=False, contig="1",
                up_shift=1000, dw_shift=1000, merged=True, interval_types="gene",
                target_genes="ENSG00000187634"):
    """Make input matrix"""
    var_header = var_hand.header
    # var_sample = var_hand.samples

    gene_matrix = {}
    for gtf_record in inter_hand.fetch("1"):
        # Important fields
        # attrbutes = gtf_record.attributes

        feature = gtf_record.feature
        if not isinstance(interval_type, list):
            interval_types = [interval_types]

        if feature not in interval_types:
            continue
        
        if not isinstance(target_genes, list):
            target_genes = [target_genes]

        gene_id = gtf_record.gene_id
        if gene_id not in target_genes:
            continue

        contig = gtf_record.contig
        strand = gtf_record.strand
        start = gtf_record.start
        end = gtf_record.end

        iv_start, iv_stop = start, end
        iv_vars = var_hand.fetch(contig=contig, start=iv_start, stop=iv_stop, reopen=True)  # 0-based indexing
        iv_vars_hash = _make_variant_dict(iv_vars)
        print(iv_vars_hash)

        # if iv_vars_hash:
            # iv_seq = seq_hand.fetch(reference=contig, start=iv_start, end=iv_stop)  # 0-based indexing

            # iv_pileups = aln_hand.pileup(contig=contig, start=iv_start, stop=iv_stop, min_base_quality=20)  # 0-based
        
            # for _column in iv_pileups:
            #     sequence_per_locus = _column.get_query_sequences()
            #     ref_pos = _column.reference_pos
            #     ref_id = _column.reference_id

            #     if (str(contig), iv_start) in iv_vars_hash:
            #         print("Ref pos: {}:{}. Also in VCF: {}".format(ref_id, ref_pos, sequence_per_locus))
            # print("CHROM: {}, position: {}", ref_id, ref_pos)

        up_matrix = None
        if up_shift and feature == "gene":
            up_start = (int(start) - 1) - up_shift
            up_end = int(start) - 1
            up_vars = var_hand.fetch(contig=contig, start=up_start, stop=up_end, reopen=True)
            up_vars_hash = _make_variant_dict(up_vars)

            if up_vars_hash:
                up_seq = seq_hand.fetch(reference=contig, start=up_start, end=up_end)
                up_matrix = matrix_factory(up_seq, up_vars_hash, contig, up_start)
                print(up_matrix)
        
        dw_matrix = None
        if dw_shift and feature == "gene":
            dw_start = int(end) + 1
            dw_end = (int(end) + 1) + dw_shift
            dw_vars = var_hand.fetch(contig=contig, start=dw_start, stop=dw_end, reopen=True)

            dw_vars_hash = _make_variant_dict(dw_vars)
            if dw_vars_hash:
                dw_seq = seq_hand.fetch(reference=contig, start=dw_start, end=dw_end)
                dw_matrix = matrix_factory(dw_seq, dw_vars_hash, contig, dw_start)
                print(dw_matrix)

        # TODO: should yield a merged metrix of upstream and downstream
        if strand == '-':
            gene_matrix[gene_id] = {"upstream": dw_matrix, "dwstream": up_matrix}
        else:
            gene_matrix[gene_id] = {"upstream": up_matrix, "dwstream": dw_matrix}

    return gene_matrix

def main(): 
    parser = get_args()
    args = parser.parse_args()

    reference = args.reference
    variants = args.variants
    annotations = args.annotations
    alignments = args.alignments

    sequence_hand = ps.FastaFile(reference)
    variant_hand = ps.VariantFile(variants, duplicate_filehandle=True)
    interval_hand = ps.TabixFile(annotations, parser=ps.asGTF())
    alignment_hand = ps.AlignmentFile(alignments, "rb", duplicate_filehandle=True)

    matrix_pool = make_matrix(interval_hand, sequence_hand, variant_hand, alignment_hand)

    sequence_hand.close()
    variant_hand.close()
    alignment_hand.close()
    interval_hand.close()

if __name__ == "__main__":
    main()


# Sequence is 0-based
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