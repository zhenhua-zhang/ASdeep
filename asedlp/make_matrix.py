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

# With direction
NT2AMB = {
    "AA": "A", "CC": "C", "GG": "G", "TT": "T", "NN": "N",
    "AC": "K", "AG": "Y", "AT": "W", "CG": "S", "CT": "R", "GT": "M",
    "CA": "k", "GA": "y", "TA": "w", "GC": "s", "TC": "r", "TG": "m",
}

# With direction
AMB2NT = {
    "A": "AA", "C": "CC", "G": "GG", "T": "TT", "N": "NN",
    "K": "AC", "Y": "AG", "W": "AT", "S": "CG", "R": "CT", "M": "GT",
    "k": "CA", "y": "GA", "w": "TA", "s": "GC", "r": "TC", "m": "TG",
}

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


def _make_variant_dict(variants):
    # Using all filter: PASS / Inaccessibl. FIXME: could be duplicated positionse
    record_dict = {(var.chrom, var.pos): var for var in variants if var.is_snp()}
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


class OrfSeqFactory(object):
    """A factory for ORF sequence.

    TODO:
        1. Not able to deal with meta-exons.
        2. Which transcript should be considered?
    """
    def __init__(self, itv_hd=None, var_hd=None, aln_hd=None):
        self.itv_hd = itv_hd
        self.var_hd = var_hd
        self.aln_hd = aln_hd

    def count(self):
        """Count reads mapped to each haplotype for given genomic interval.
        """
        contig, start, end = self.itv_hd.contig, self.itv_hd.start, self.itv_hd.end
        itv_vars = self.var_hd.fetch(contig=contig, start=start, stop=end, reopen=True)
        itv_vars_hash = _make_variant_dict(itv_vars)

        if itv_vars_hash:
            gene_id = self.itv_hd.gene_id
            pileup_seg = self.aln_hd.pileup(contig=contig, start=start, stop=end)
            for pileup_col in pileup_seg:
                ref_pos = pileup_col.reference_pos
                # pileups = pileup_col.pileups
                qry_bases = pileup_col.get_query_sequences(mark_matches=True, mark_ends=True, add_indels=True)
        else:
            print("No heterogeous loci in: {}:{}-{}".format(contig, start, end))

    def binomial_test(self):
        """Do Binomial test on given data.
        """

    def beta_binomial_test(self):
        """Do Beta-binomial test on given data.
        """

    def ase_test(self, mthd="bb"):
        """Determine whether there is an ASE effects.
        """
        return True


class NonOrfSeqFactory(object):
    """A factory for non-ORF sequence.
    """
    def __init__(self, itv_hd=None, seq_hd=None, var_hd=None, aln_hd=None):
        self.itv_hd = itv_hd
        self.seq_hd = seq_hd
        self.var_hd = var_hd
        self.aln_hd = aln_hd

    def seq2matrix(self, seq, var_hash, chrom, shift, target_samples):
        """Make a sequence matrix of upstream of a gene.
        TODO:
            1. Should `shift` be determined by strand?
        """
        if not isinstance(var_hash, dict):
            raise TypeError("Require dict")

        _seq_init_len = len(seq)
        target_sample_allele_vec = {}
        for _each_sample in target_samples: # Using one thread for each sample.
            for (_chrom, _pos), _vcf_rec in var_hash.items(): # multiple threading could be useful
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
                        .get(_each_sample, {"GT": (0, 0)}) \
                        .get("GT", (0, 0))
                _sub_base = NT2AMB[_alleles[_phase_0] + _alleles[_phase_1]]
                seq = seq[:_pos] + _sub_base + seq[_pos+1:]

            _seq_final_len = len(seq)

            if _seq_init_len != _seq_final_len:
                continue

            _allele_vec = []
            for _base in seq:
                _allele_0, _allele_1 = AMB2NT.get(_base, "NN")
                _allele_vec.append(NT2VEC.get(_allele_0, DNTVEC) + NT2VEC.get(_allele_1, DNTVEC))

            target_sample_allele_vec[_each_sample] = _allele_vec

        return target_sample_allele_vec
    
    def enc_nc_reg(self, contig, start, end, strand, dw_shift, up_shift, target_samples="gonl-100a"):
        """Encode regutlation sequence into matrix.
        """
        if isinstance(target_samples, str):
            target_samples = [target_samples]

        up_matrix = None
        if up_shift:
            up_start = (int(start) - 1) - up_shift
            up_end = int(start) - 1
            up_vars = self.var_hd.fetch(contig=contig, start=up_start, stop=up_end, reopen=True)
            up_vars_hash = _make_variant_dict(up_vars)

            if up_vars_hash:  # If there's heterogeous loci in upstream
                up_seq = self.seq_hd.fetch(reference=contig, start=up_start, end=up_end)
                up_matrix = self.seq2matrix(up_seq, up_vars_hash, contig, up_start, target_samples)
        
        dw_matrix = None
        if dw_shift:
            dw_start = int(end) + 1
            dw_end = (int(end) + 1) + dw_shift
            dw_vars = self.var_hd.fetch(contig=contig, start=dw_start, stop=dw_end, reopen=True)

            dw_vars_hash = _make_variant_dict(dw_vars)
            if dw_vars_hash:  # If there's heterogeous loci in downstream
                dw_seq = self.seq_hd.fetch(reference=contig, start=dw_start, end=dw_end)
                dw_matrix = self.seq2matrix(dw_seq, dw_vars_hash, contig, dw_start, target_samples)

        if strand == '-':
            return {"upstream": dw_matrix, "dwstream": up_matrix}
        else:
            return {"upstream": up_matrix, "dwstream": dw_matrix}


class ASEFactory:
    """A factory to work on ASE related files.
    """
    def __init__(self, itv_hand=None, seq_hand=None, var_hand=None, aln_hand=None):
        self.itv_hd = itv_hand
        self.seq_hd = seq_hand
        self.var_hd = var_hand
        self.aln_hd = aln_hand

    def additvhd(self, itv_file):
        self.itv_hd = ps.TabixFile(itv_file, parser=ps.asGTF())
        return self

    def addseqhd(self, seq_file):
        self.seq_hd = ps.FastaFile(seq_file)
        return self

    def addvarhd(self, var_file):
        self.var_hd = ps.VariantFile(var_file, duplicate_filehandle=True)
        return self

    def addalnhd(self, aln_file):
        self.aln_hd = ps.AlignmentFile(aln_file, "rb", duplicate_filehandle=True)
        return self

    def shutdown(self):
        if self.itv_hd:
            self.itv_hd.close()

        if self.seq_hd:
            self.seq_hd.close()

        if self.var_hd:
            self.var_hd.close()

        if self.aln_hd:
            self.aln_hd.close()

    def transform(self, contig="1", up_shift=100, dw_shift=100, target_itvl="ENSG00000187634"):
        """Encode given genomic region.

        """
        if not isinstance(target_itvl, list):
            target_itvl = [target_itvl]

        for itvrc in self.itv_hd.fetch(contig):  # `multiple_iterators` helps create non-consuming iterators
            if not itvrc:
                continue

            gene_id = itvrc.gene_id
            if gene_id not in target_itvl:
                continue

            if len(target_itvl) == 0:
                break

            target_itvl.remove(gene_id)

            feature, start, end, strand = itvrc.feature, itvrc.start, itvrc.end, itvrc.strand
            if feature == "gene": # to parse upstream- and downstream-sequence into a sparse matrix
                nosf = NonOrfSeqFactory(
                    aln_hd=self.aln_hd, itv_hd=self.itv_hd,
                    seq_hd=self.seq_hd, var_hd=self.var_hd
                )
                _reg_matrix = nosf.enc_nc_reg(contig, start, end, strand, dw_shift, up_shift)
            elif feature == "exon":  # Make haplotype read count matrix for each exon
                osf = OrfSeqFactory(itvrc, self.var_hd, self.aln_hd)
                _ase_effect = osf.ase_test()
            else:
                print("Skip feature: {} ...".format(feature))
            print(_ase_effect, _reg_matrix)


def main(): 
    """Main function.
    """
    args = get_args().parse_args()

    variants = args.variants
    reference = args.reference
    alignments = args.alignments
    annotations = args.annotations

    factory = ASEFactory()
    factory.addalnhd(alignments).addseqhd(reference).additvhd(annotations).addvarhd(variants)

    _matrix_pool = factory.transform()
    print(_matrix_pool)

    factory.shutdown()


if __name__ == "__main__":
    main()

