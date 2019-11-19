#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Make matrix"""
import os
import sys

from argparse import ArgumentParser

import pysam as ps

DNTVEC = [0, 0, 0, 0]  # Base not in ACGTN
NT2VEC = {"A": [1, 0, 0, 0], "C": [0, 1, 0, 0], "G": [0, 0, 1, 0], "T": [0, 0, 0, 1], "N": [0, 0, 0, 0]}
VEC2NT = {(1, 0, 0, 0): "A", (0, 1, 0, 0): "C", (0, 0, 1, 0): "G", (0, 0, 0, 1): "T", (0, 0, 0, 0): "N"}

NT2AMB = {
    "AA": "A", "CC": "C", "GG": "G", "TT": "T", "NN": "N",
    "AC": "K", "AG": "Y", "AT": "W", "CG": "S", "CT": "R", "GT": "M",
    "CA": "k", "GA": "y", "TA": "w", "GC": "s", "TC": "r", "TG": "m",
}

AMB2NT = {
    "A": "AA", "C": "CC", "G": "GG", "T": "TT", "N": "NN",
    "K": "AC", "Y": "AG", "W": "AT", "S": "CG", "R": "CT", "M": "GT",
    "k": "CA", "y": "GA", "w": "TA", "s": "GC", "r": "TC", "m": "TG",
}

def get_args():
    """Get CLI arguments for current script"""
    parser = ArgumentParser()
    _group = parser.add_argument_group("Input")
    _group.add_argument(
        "-r", "--reference", type=str, dest="reference", action="store",
        required=True, help="Reference genome (FASTA). Required"
    )
    _group.add_argument(
        "-v", "--variants", type=str, dest="variants", action="store",
        required=True, help="Vairant file (VCF). Required"
    )
    _group.add_argument(
        "-a", "--annotations", type=str, dest="annotations", action="store",
        required=True, help="Annotation file in (GFF/GTF). Required"
    )
    _group.add_argument(
        "-s", "--alignments", type=str, dest="alignments", action="store",
        required=True, help="Sequence alignment file (SAM/BAM). Required"
    )
    return parser


class UnknownStrandError(Exception):
    pass


class LengthError(Exception):
    pass


class NonDictError(Exception):
    pass


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


class Binomial(object):
    pass


class BetaBinomial(object):
    pass


class HapReadsCounter(object):
    """Count reads mapped to haplotypes"""

    def __init__(self, itv_recd, var_hand, aln_hand):
        self.itv_recd = itv_recd
        self.var_hand = var_hand
        self.aln_hand = aln_hand

    def count(self):
        """Count reads mapped to each haplotype for given genomic interval"""
        # TODO: perhaps need exon_id as well, e.g. smooth meta-exon
        end = self.itv_recd.end
        start = self.itv_recd.start
        contig = self.itv_recd.contig

        itv_vars = self.var_hand.fetch(contig=contig, start=start, stop=end, reopen=True)
        itv_vars_hash = _make_variant_dict(itv_vars)

        if itv_vars_hash:
            gene_id = self.itv_recd.gene_id
            pileup_seg = self.aln_hand.pileup(contig=contig, start=start, stop=end)
            for pileup_col in pileup_seg:
                ref_pos = pileup_col.reference_pos
                # pileups = pileup_col.pileups
                qry_bases = pileup_col.get_query_sequences(mark_matches=True, mark_ends=True, add_indels=True)
                print(gene_id)
                print(ref_pos)
                print(qry_bases)
        else:
            print("No heterogeous loci in: {}:{}-{}".format(contig, start, end))

    def _beta_binom_test(self):
        """Do Beta-binomial test on given data"""

    def _binom_test(self):
        """Do Binomial test on given data"""

    def check_ase(self, mthd="bb"):
        """Determine whether there is an ASE effects"""

    def concat(self, counter):
        """Concatenate current counter with another"""


class SequenceMatrixFactory(object):
    """Sequence matrix of regulation region"""

    def __init__(self):
        self.itv_hand = None
        self.seq_hand = None
        self.var_hand = None
        self.aln_hand = None

    def new(self):
        """Create a new sequence matrix from given genomic interval, genomic variants"""
        pass

    def into_narray(self):
        """Transform the matrix into a numpy"""

    def get_region_string(self):
        pass

    def get_current_coord(self):
        """Return the coordination of current sequence matrix"""

    def add_itv_handle(self, itv_file):
        self.itv_hand = ps.TabixFile(itv_file, parser=ps.asGTF())
        return self

    def add_seq_handle(self, seq_file):
        self.seq_hand = ps.FastaFile(seq_file)
        return self

    def add_var_handle(self, var_file):
        self.var_hand = ps.VariantFile(var_file, duplicate_filehandle=True)
        return self

    def add_aln_handle(self, aln_file):
        self.aln_hand = ps.AlignmentFile(aln_file, "rb", duplicate_filehandle=True)
        return self

    def close_all_handle(self):
        if self.itv_hand:
            self.itv_hand.close()

        if self.seq_hand:
            self.seq_hand.close()

        if self.var_hand:
            self.var_hand.close()

        if self.aln_hand:
            self.aln_hand.close()

    def _generate_sequence_matrix(self, seq, var_hash, chrom, shift,
            target_samples=["gonl-100a"]):
        """Make a sequence matrix of upstream of a gene"""
        # FIXME: `shift` should be determined by strand???
        if not isinstance(var_hash, dict):
            raise NonDictError  # TODO: a concrete sub-class of TypeError

        if isinstance(target_samples, str):
            target_samples = [target_samples]

        _seq_init_len = len(seq)
        target_sample_allele_vec = {}
        for _each_sample in target_samples:
            for (_chrom, _pos), _vcf_rec in var_hash.items():
                _pos = _pos - shift - 1  # TODO: Not sure about the coordination

                _alleles = _vcf_rec.alleles  # tuple of reference allele followed by alt alleles
                _allele_a, _allele_b = _alleles
                if len(_allele_a) != 1 or len(_allele_b) != 1:  # indels will be skipped
                    continue

                _ref_allele = seq[_pos]  # Reference allele base in sequence
                if _ref_allele != _allele_a:
                    continue

                # If there is not a GT, will use the reference allele
                _phase_0, _phase_1 = _vcf_rec.samples \
                    .get(_each_sample, None) \
                    .get("GT", (0, 0))
                _sub_base = NT2AMB[_alleles[_phase_0] + _alleles[_phase_1]]
                seq = seq[:_pos] + _sub_base + seq[_pos+1:]

            _seq_final_len = len(seq)

            if _seq_init_len != _seq_final_len:
                continue

            _allele_vec = []
            for _base in seq:
                _allele_0, _allele_1 = AMB2NT.get(_base, "NN")
                _allele_vec.append(
                    NT2VEC.get(_allele_0, DNTVEC) +
                    NT2VEC.get(_allele_1, DNTVEC)
                )

            target_sample_allele_vec[_each_sample] = _allele_vec

        # TODO: perhaps a Pandas DataFrame will make life easier
        return target_sample_allele_vec

    def _count_haplotype_read(self, itv_recd):
        """Count reads mapped to each haplotype """
        counter = HapReadsCounter(itv_recd, self.var_hand, self.aln_hand)
        counter.count()

    def make_matrix(self, with_orf=False, contig="1", up_shift=100,
            dw_shift=100, target_genes="ENSG00000187634"):
        """Make input matrix"""
        if not isinstance(target_genes, list):
            target_genes = [target_genes]

        gene_matrix = {}
        for itv_recd in self.itv_hand.fetch(contig):  # `multiple_iterators` helps create non-consuming iterators
            if not itv_recd:
                continue

            gene_id = itv_recd.gene_id
            if gene_id not in target_genes:
                continue

            feature = itv_recd.feature
            strand, start, end = itv_recd.strand, itv_recd.start, itv_recd.end

            if feature == "gene": # to parse upstream- and downstream-sequence into a parse matrix
                up_matrix = None
                if up_shift and feature == "gene":
                    up_start = (int(start) - 1) - up_shift
                    up_end = int(start) - 1
                    up_vars = self.var_hand.fetch(contig=contig, start=up_start, stop=up_end, reopen=True)
                    up_vars_hash = _make_variant_dict(up_vars)

                    if up_vars_hash:  # If there's heterogeous loci in upstream
                        up_seq = self.seq_hand.fetch(reference=contig, start=up_start, end=up_end)
                        up_matrix = self._generate_sequence_matrix(up_seq, up_vars_hash, contig, up_start)
                
                dw_matrix = None
                if dw_shift and feature == "gene":
                    dw_start = int(end) + 1
                    dw_end = (int(end) + 1) + dw_shift
                    dw_vars = self.var_hand.fetch(contig=contig, start=dw_start, stop=dw_end, reopen=True)

                    dw_vars_hash = _make_variant_dict(dw_vars)
                    if dw_vars_hash:  # If there's heterogeous loci in downstream
                        dw_seq = self.seq_hand.fetch(reference=contig, start=dw_start, end=dw_end)
                        dw_matrix = self._generate_sequence_matrix(dw_seq, dw_vars_hash, contig, dw_start)

                # TODO: should yield a merged metrix of upstream and downstream
                if strand == '-':
                    gene_matrix[gene_id] = {"upstream": dw_matrix, "dwstream": up_matrix}
                else:
                    gene_matrix[gene_id] = {"upstream": up_matrix, "dwstream": dw_matrix}
                    pass
            elif feature == "exon":  # Make haplotype read count matrix for each exon
                self._count_haplotype_read(itv_recd)
            else:
                pass

        return self

    # def het_reads_counter(self, itv_rec, var_hash, aln_hand):
    #     """Quantify ASE effects by heterozygous loci in exon regions for a gene
    #     @param itv_rec Genomic intervl record
    #     @param var_hash Variants dictionary
    #     @param aln_hand Sequence alignment file handle
    #     """
    #     pass


def main(): 
    parser = get_args()
    args = parser.parse_args()

    reference = args.reference
    variants = args.variants
    annotations = args.annotations
    alignments = args.alignments

    sequence_matrix = SequenceMatrixFactory()
    sequence_matrix \
        .add_aln_handle(alignments) \
        .add_seq_handle(reference) \
        .add_itv_handle(annotations) \
        .add_var_handle(variants)

    matrix_pool = sequence_matrix.make_matrix()

    sequence_matrix.close_all_handle()


if __name__ == "__main__":
    main()


class SequenceMatrixPool(object):
    """A container including multiple SequenceMatrix"""

    def __init__(self):
        pass

    def add_matrix(self, sequence_matrix):
        """Add an new matrix into current SequenceMatrixPool"""
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
