#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import gzip
import pprint
import argparse

from collections import UserDict

from scipy.stats import binom, chi2, betabinom
from scipy.optimize import minimize

import tables
import gffutils
import numpy as np


def get_args():
    parser = argparse.ArgumentParser(
        prog="ASEFactory", description="A factory to produce ASE effects.")

    parser.add_argument(
        "-H", "--haplotypes", action="store", default="./haps.h5",
        type=str, dest="haplotypes",
        help="Path to HDF5 file to read phased haplotypes from."
    )

    parser.add_argument(
        "-i", "--snp-index", action="store", default="./snp_idex.h5",
        type=str, dest="snp_index",
        help="Path to HDF5 file to read SNP index from."
    )

    parser.add_argument(
        "-s", "--snp-tab", action="store", default="./snp_tab.h5",
        type=str, dest="snp_tab",
        help="Path to HDF5 file to read SNP information from."
    )

    parser.add_argument(
        "-S", "--sequence-tab", action="store", default="./seq_tab.h5",
        type=str, dest="seq_tab",
        help="Path to HDF5 file to read reference sequence from."
    )

    parser.add_argument(
        "-r", "--ref-read-counts", action="store", default="./ref_count_tab.h5",
        type=str, dest="ref_tab",
        help="Path to HDF5 file to read reference reads counts from."
    )

    parser.add_argument(
        "-a", "--alt-read-counts", action="store", default="./alt_count_tab.h5",
        type=str, dest="alt_tab",
        help="Path to HDF5 file to read alternative reads counts from."
    )

    parser.add_argument(
        "-g", "--genome-annot", action="store", default="./genome.gff.gz",
        type=str, dest="genome_annot",
        help="Path to GFF / GTF file to read gene structure information from."
    )

    return parser.parse_args()


class ReadCountPool(UserDict):
    """A class for Read counts.
    """
    def __init__(self):
        super(ReadCountPool, self).__init__()
        self.data = {}

    def __setitem__(self, key, value):
        if key in self.data:
            self.data[key] = value
        else:
            self.data[key] = value

    def __getitem__(self, key, default=None):
        if key in self.data:
            return self.data[key]
        else:
            return default

    def __len__(self):
        return len(self.data)

    def _p_add_counts(self, id, counts, type="oth"):
        """Add counts.
        """
        type_pool = ["ref", "alt"]
        if id not in self.data:
            self.data[id] = [[], []]

        self.data[id][type_pool.index(type)] = counts

        len_pool = [len(x) for x in self.data[id]]
        if 0 not in len_pool:
            assert len_pool[0] == len_pool[1]

    def add_ref_counts(self, id, counts):
        self._p_add_counts(id, counts, type="ref")

    def add_alt_counts(self, id, counts):
        self._p_add_counts(id, counts, type="alt")

    def unpack(self, _mthd=None):
        """Unpack the object.
        """
        # TODO: Add a _mthd parameter to give more control of the unpacked vale.
        ref_read_pool, alt_read_pool = [], []
        for _key, _value in self.data.items():
            ref_read_pool.extend(_value[0])
            alt_read_pool.extend(_value[1])
        return ref_read_pool, alt_read_pool


class ASEFactory:
    """A class to produce ASE effects and matrix of regulatory sequence"""
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

    NNTVEC = [0, 0, 0, 0]  # Base not in ACGTN
    NT2VEC = {"A": [1, 0, 0, 0], "C": [0, 1, 0, 0], "G": [0, 0, 1, 0], "T": [0, 0, 0, 1]}
    VEC2NT = {(1, 0, 0, 0): "A", (0, 1, 0, 0): "C", (0, 0, 1, 0): "G", (0, 0, 0, 1): "T"}

    _sample_id_to_idx = None
    exon_pool = None
    hap_tab, snp_idx, snp_tab, seq_tab, ref_tab, alt_tab, ant_sql = [None] * 7

    def __init__(self, args):
        super(ASEFactory, self).__init__()
        self.args = args
        self.ase = None
        self.mrna = None
        self.ntmtrx_pool = None

        self._p_check_args()

    def _p_check_args(self):
        _args = self.args

    def _p_sample_id_to_idx(self, chrom, sample_id="gonl-100a"):
        if not isinstance(sample_id, list):
            sample_id = [sample_id]

        if self._sample_id_to_idx is None:
            hap_samples = self.hap_tab.get_node("/samples_{}".format(chrom))
            self._sample_id_to_idx = {
                _sample[0].decode(): _idx
                for _idx, _sample in enumerate(hap_samples)
            }

        return [self._sample_id_to_idx.get(_id) for _id in sample_id]

    def _p_get_nodes(self, chrom, requested_nodes=("haplotypes", "read_counts", "SNPs", "sequences")):
        node_pool = []
        for node_name in requested_nodes:
            if "haplotypes" == node_name:
                node_pool.append(self.hap_tab.get_node("/{}".format(chrom)))
                node_pool.append(self.hap_tab.get_node("/phase_{}".format(chrom)))
            elif "read_counts" == node_name:
                node_pool.append(self.alt_tab.get_node("/{}".format(chrom)))
                node_pool.append(self.ref_tab.get_node("/{}".format(chrom)))
            elif "SNPs" == node_name:
                node_pool.append(self.snp_tab.get_node("/{}".format(chrom)))
                node_pool.append(self.snp_idx.get_node("/{}".format(chrom)))
            elif "sequences" == node_name:
                node_pool.append(self.seq_tab.get_node("/{}".format(chrom)))
            else:
                print("Unsupported type of node: {}".format(node_name))

        return tuple(node_pool)

    def _p_nt2vc(self, code):
        if isinstance(code, (np.int8, int)):
            code = chr(code)
        elif isinstance(code, bytes):
            code = code.decode()
        return self.NT2VEC.get(code, self.NNTVEC)

    def _p_encode_nt(self, a1_code, a2_code):
        """Encode allele pairs into matrix"""
        return self._p_nt2vc(a1_code) + self._p_nt2vc(a2_code)

    def _p_try_open(self, file_path):
        try:
            return tables.open_file(file_path)
        except FileExistsError as feerr:
            print(feerr, file=sys.stderr)
        except PermissionError as pmerr:
            print(pmerr, file=sys.stderr)
        except IOError as ioerr:
            print(ioerr, file=sys.stderr)
        except Exception as exp:
            print(exp, file=sys.stderr)

    def _p_try_close(self, tables_obj):
        try:
            tables_obj.close()
        except IOError as ioerr:
            print(ioerr, file=sys.stderr)
        except Exception as exp:
            print(exp, file=sys.stderr)

    def init(self, ant_db_name=":memory:"):
        args = self.args
        self.hap_tab = self._p_try_open(args.haplotypes)
        self.snp_idx = self._p_try_open(args.snp_index)
        self.snp_tab = self._p_try_open(args.snp_tab)
        self.seq_tab = self._p_try_open(args.seq_tab)
        self.ref_tab = self._p_try_open(args.ref_tab)
        self.alt_tab = self._p_try_open(args.alt_tab)

        if not os.path.exists(ant_db_name):
            gffutils.create_db(args.genome_annot, ant_db_name,
                               disable_infer_transcripts=True,
                               disable_infer_genes=True)

        self.ant_sql = gffutils.FeatureDB(ant_db_name)
        return self

    def shutdown(self):
        # TODO: if the attribute is None it shouldn't have a close() method
        self.hap_tab.close()
        self.snp_idx.close()
        self.snp_tab.close()
        self.ref_tab.close()
        self.alt_tab.close()
        self.seq_tab.close()

    def gen_gnm_region(self, gene_id, featuretype="transcript", which_mrna=None):
        """Fetches exon regions for given gene ID.

        NOTE:
            A more careful calibration is required to choose mRNA.
        Args:
            gene_id (str):
            featuretype (str, optional, None):
            which_mrna (str, optional, None): The message RNA id
        """
        mrna_pool = self.ant_sql.children(gene_id, featuretype=featuretype)

        if which_mrna is None:
            mrna = next(mrna_pool, False)
        else:
            raise NotImplementedError("Not implemented yet")

        if mrna:
            parent_id = mrna.id
        else:
            parent_id = gene_id

        exon_pool = self.ant_sql.children(parent_id, featuretype="exon")

        self.mrna = mrna
        self.exon_pool = exon_pool

        return self

    def gen_seq_mtrx(self, seq_itvl=None, sample_id="gonl-100a", with_dw=False, shift=1e2):
        """Generate a chain of amb for a sequence.

        TODO:
            1. Allow different size of shrinkage for up- and down-stream.
        """
        if not isinstance(sample_id, list):
            sample_id = [sample_id]

        if seq_itvl is None:
            seq_itvl = self.mrna

        chrom = seq_itvl.seqid

        sample_idx = self._p_sample_id_to_idx(chrom, sample_id)
        (seq_code_pool, snp_indx_pool, snp_code_pool, hap_code_pool,
         hap_phase_pool) = self._p_get_nodes(chrom, ("sequences", "SNPs", "haplotypes"))

        anchor_pool = [[int(seq_itvl.start - shift), seq_itvl.start]]
        if with_dw:  # FIXME: Not sure if with the end base
            anchor_pool = anchor_pool.append([seq_itvl.stop + 1, seq_itvl.stop + 1 + shift])

        ntmtrx_pool = {0: {}, 1: {}}
        for anchor_idx, anchor in enumerate(anchor_pool):  # No more than twice
            seq_code = seq_code_pool[anchor[0]: anchor[1]]
            snp_code = snp_code_pool[anchor[0]: anchor[1]]
            for _id, _idx in zip(sample_id, sample_idx):  # No more than # of sample_id
                ntmtrx = []
                for a1_code, a2_code in list(zip(seq_code, snp_code)):
                    if a2_code != -1:
                        if hap_phase_pool[a2_code, sample_idx] == 1:  # Is phased
                            snp_info = snp_indx_pool[a2_code]
                            a1_idx, a2_idx = hap_code_pool[a2_code, _idx:_idx + 2]
                            a1_code, a2_code = snp_info[a1_idx + 2], snp_info[a2_idx + 2]
                    ntmtrx.append(self._p_encode_nt(a1_code, a1_code))
                if _id in ntmtrx_pool[anchor_idx]:
                    raise KeyError(_id + " shouldn't in the ntmtrx_pool")
                ntmtrx_pool[anchor_idx][_id] = ntmtrx
        self.ntmtrx_pool = ntmtrx_pool

        return self

    def gen_ase_effect(self, seq_itvl=None, exon_pool=None, mthd="bn", meta_exon=False):
        """Generate allele-specific expression effects from read counts.
        """
        if seq_itvl is None:
            seq_itvl = self.mrna

        chrom = seq_itvl.seqid
        ref_read_pool, alt_read_pool = self._p_get_nodes(chrom, ("read_counts", ))

        if exon_pool is None:
            exon_pool = self.exon_pool

        rcp = ReadCountPool()
        for exon_itvl in exon_pool:
            exon_start, exon_stop = exon_itvl.start, exon_itvl.stop
            id = (exon_itvl.seqid, exon_start, exon_stop) 
            ref_counts = ref_read_pool[exon_start: exon_stop]
            alt_counts = alt_read_pool[exon_start: exon_stop]

            rcp.add_ref_counts(id, ref_counts)
            rcp.add_alt_counts(id, alt_counts)

        aefact = ASEeffectFactory(rcp)
        if mthd == "bb":
            self.ase = aefact.bbtest()
        elif mthd == "bn":
            self.ase = aefact.bntest()

        return self
    
    def report(self):
        print("Region id: ", self.mrna.id)
        print("Likelihood: ", self.ase[0])
        print("P-value: ", self.ase[1])
        print("Sequence matrix: ", self.ntmtrx_pool)
        return self


class ASEeffectFactory:
    """A factory for ORF sequence.

    TODO:
        1. Not able to deal with meta-exons.
        2. Which transcript to be considered.
    """
    def __init__(self, read_counts):
        self.rrc, self.arc = self._p_parse_read_counts(read_counts)
        self.opt_results = {}

    def _p_parse_read_counts(self, read_counts):
        if isinstance(read_counts, (tuple, list)) and len(read_counts) == 2:
            return read_counts
        elif isinstance(read_counts, ReadCountPool):
            return read_counts.unpack()
        else:
            raise TypeError("read_counts parameter should be either tuple, list, or dict")

    def _p_bnllh(self, param, k_vec, n_vec):
        """The likelihood function of Binomial distribution.
        """
        lh_vec = [binom.logpmf(k, n, param) for k, n in zip(k_vec, n_vec)]
        return -sum(lh_vec)

    def bntest(self):
        """Do Binomial test on given data.
        """
        x0 = 0.5
        k_vec = self._p_get_k_vec()
        n_vec = self._p_get_n_vec()

        opt_results = minimize(
            self._p_bnllh, x0, args=(k_vec, n_vec),
            method="nelder-mead", options={"maxfev": 1e3, "ftol": 1e-8}
        )
        self.opt_results["binomial"] = opt_results

        _est_t = opt_results["x"]
        estll = opt_results["fun"]

        hyp_t = 0.5
        hypll = self._p_bnllh(hyp_t, k_vec, n_vec)

        llr = -2 * (estll - hypll)
        p = chi2.sf(llr, 1)

        return llr, p

    def _p_bbllh(self, params, k_vec, n_vec):
        """The likelihood function of Beta-Binomial distribution.
        """
        if len(params) != 2 or not isinstance(params, list):
            raise TypeError("The params should be a list with two elements.")

        bb_a, bb_b = params
        lh_vec = [betabinom.logpmf(k, n, bb_a, bb_b) for k, n in zip(k_vec, n_vec)]
        return -sum(lh_vec)

    def bbtest(self):
        """Do Beta-binomial test on given data.
        """
        x0 = [0.5, 0.5]
        k_vec = self._p_get_k_vec()
        n_vec = self._p_get_n_vec()

        opt_results = minimize(
            self._p_bbllh, x0, args=(k_vec, n_vec),
            method="nelder-mead", options={"maxfev": 1e3, "ftol": 1e-8}
        )
        self.opt_results["beta_binomial"] = opt_results

        _est_a, _est_b = opt_results["x"]
        estll = opt_results["fun"]

        hyp_a = hyp_b = opt_results["x"].mean()
        hypll = self._p_bbllh((hyp_a, hyp_b), k_vec, n_vec)

        llr = - 2 * (estll - hypll) # The likelihood value is timed by -1
        p = chi2.sf(llr, 1)

        return llr, p

    def _p_get_k_vec(self):
        return self.rrc

    def _p_get_n_vec(self):
        return list(map(sum, zip(self.rrc, self.arc)))

    def chi_square_test(self):
        """Using a Chi-square test on given data.
        """
        return NotImplemented


def main():
    args = get_args()
    factory = ASEFactory(args)
    factory.init("test.db") \
            .gen_gnm_region(gene_id="ENSG00000142657") \
            .gen_seq_mtrx() \
            .gen_ase_effect() \
            .report() \
            .shutdown()


if __name__ == '__main__':
    main()

