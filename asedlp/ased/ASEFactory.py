"""ASE factory."""

import os
import sys
import math
from collections import UserDict

import tables
import gffutils
import numpy as np
from scipy.optimize import minimize
from scipy.stats import binom, chi2, betabinom

from .zutils import UserDict
from .zutils import cmp


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

    def _p_add(self, itvl_id: str, content: list, add_type="a1"):
        """Add counts.
        """
        type_pool = ["a1", "a2", "info"]
        if itvl_id not in self.data:
            self.data[itvl_id] = [[], [], []]

        self.data[itvl_id][type_pool.index(add_type)] = content

    def add_a1_rc(self, itvl_id, counts):
        self._p_add(itvl_id, counts, add_type="a1")

    def add_a2_rc(self, itvl_id, counts):
        self._p_add(itvl_id, counts, add_type="a2")

    def add_loci_info(self, itvl_id, info):
        self._p_add(itvl_id, info, add_type="info")

    def unpack(self):
        """Unpack the object.
        """
        a1_rc_pool, a2_rc_pool, info_pool = [], [], []
        for _key, _value in self.data.items():
            a1_rc_pool.extend(_value[0])
            a2_rc_pool.extend(_value[1])
            info_pool.extend(_value[2])
        return a1_rc_pool, a2_rc_pool, info_pool


class ASEeffectFactory:
    """A factory for ORF sequence.

    TODO:
        1. Not able to deal with meta-exons.
        2. Which transcript to be considered.
    """

    def __init__(self, rc):
        self.rc = self._p_parse_rc(rc)
        self.optim_results = {}

    @staticmethod
    def _p_parse_rc(rc):
        if isinstance(rc, (tuple, list)) and len(rc) >= 2:
            return rc, []
        elif isinstance(rc, ReadCountPool):
            return rc.unpack()
        else:
            raise TypeError("rc parameter should be tuple, list, or ReadCountPool")

    @staticmethod
    def _p_bnllh(param, k_vec, n_vec):
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
        self.optim_results["binomial"] = opt_results

        _est_t = opt_results["x"]
        estll = opt_results["fun"]

        hyp_t = 0.5
        hypll = self._p_bnllh(hyp_t, k_vec, n_vec)

        llr = - 2 * (estll - hypll)
        p = chi2.sf(llr, 1)

        k_sum, n_sum = int(sum(k_vec)), int(sum(n_vec))
        return llr, p, cmp(2 * k_sum, n_sum)  # likelihood ratio, p-value, direction

    @staticmethod
    def _p_bbllh(params, k_vec, n_vec):
        """The likelihood function of Beta-Binomial distribution.
        """
        if len(params) != 2:
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
        self.optim_results["beta_binomial"] = opt_results

        _est_a, _est_b = opt_results["x"]
        estll = opt_results["fun"]

        hyp_a = hyp_b = opt_results["x"].mean()
        hypll = self._p_bbllh([hyp_a, hyp_b], k_vec, n_vec)

        llr = - 2 * (estll - hypll)  # The likelihood value is timed by -1
        p = chi2.sf(llr, 1)

        k_sum, n_sum = int(sum(k_vec)), int(sum(n_vec))
        return llr, p, cmp(2 * k_sum, n_sum)  # likelihood ratio, p-value, direction

    def _p_get_k_vec(self):
        return self.rc[0]

    def _p_get_n_vec(self):
        return list(map(sum, zip(self.rc[0], self.rc[1])))

    def chi_square_test(self):
        """Using a Chi-square test on given data.
        """
        return NotImplemented


class ASEFactory:
    """A class to produce ASE effects and matrix of regulatory sequence.
    """
    NT2AMB = {"AA": "A", "CC": "C", "GG": "G", "TT": "T", "NN": "N", "AC": "K", "AG": "Y", "AT": "W", "CG": "S", "CT": "R", "GT": "M", "CA": "k", "GA": "y", "TA": "w", "GC": "s", "TC": "r", "TG": "m", }
    AMB2NT = {"A": "AA", "C": "CC", "G": "GG", "T": "TT", "N": "NN", "K": "AC", "Y": "AG", "W": "AT", "S": "CG", "R": "CT", "M": "GT", "k": "CA", "y": "GA", "w": "TA", "s": "GC", "r": "TC", "m": "TG", }
    NNTVEC = [0, 0, 0, 0]  # Base not in ACGTN
    NT2VEC = {"A": [1, 0, 0, 0], "C": [0, 1, 0, 0], "G": [0, 0, 1, 0], "T": [0, 0, 0, 1]}
    VEC2NT = {(1, 0, 0, 0): "A", (0, 1, 0, 0): "C", (0, 0, 1, 0): "G", (0, 0, 0, 1): "T"}

    _SAMPLE_ID_TO_IDX = None  # Has to be None for there's a if-else branch
    _CHROM = None  # Static
    _SAMPLE_ID = None  # Static
    _gene_ids = None
    hap_tab, snp_idx, snp_tab, seq_tab, ref_tab, alt_tab, ant_sql = [None] * 7

    def __init__(self, args):
        super(ASEFactory, self).__init__()
        self.args = args
        self.mrna_pool = {}
        self.exon_pool = {}
        self.ase_pool = None
        self.ntmtrx_pool = None

    def _parse_gene_ids(self):
        self._gene_ids = self.args.gene_ids + self._parse_gene_id_file()

    def _parse_gene_id_file(self):
        if self.args.gene_id_file:
            with open(self.args.gene_id_file) as gifh:
                ids = [x.strip("\n") for x in gifh]
            return ids
        return []

    def _pr_sample_id_to_idx(self, chrom, sample_id="gonl-101b"):
        if self._SAMPLE_ID_TO_IDX is None:
            hap_samples = self.hap_tab.get_node("/samples_{}".format(chrom))
            self._SAMPLE_ID_TO_IDX = {
                _sample[0].decode(): _idx
                for _idx, _sample in enumerate(hap_samples)
            }

        return self._SAMPLE_ID_TO_IDX.get(sample_id)

    def _pr_get_nodes(self, chrom, requested_nodes=("hap", "rc", "snp", "seq")):
        node_pool = []
        for node_name in requested_nodes:
            if "hap" == node_name:
                node_pool.append(self.hap_tab.get_node("/{}".format(chrom)))
                node_pool.append(self.hap_tab.get_node("/phase_{}".format(chrom)))
            elif "rc" == node_name:
                node_pool.append(self.alt_tab.get_node("/{}".format(chrom)))
                node_pool.append(self.ref_tab.get_node("/{}".format(chrom)))
            elif "snp" == node_name:
                node_pool.append(self.snp_tab.get_node("/{}".format(chrom)))
                node_pool.append(self.snp_idx.get_node("/{}".format(chrom)))
            elif "seq" == node_name:
                node_pool.append(self.seq_tab.get_node("/{}".format(chrom)))
            else:
                print("Unsupported type of node: {}".format(node_name))

        return tuple(node_pool)

    def _pr_nt2vc(self, code):
        if isinstance(code, (np.int8, int)):
            code = chr(code)
        elif isinstance(code, bytes):
            code = code.decode()
        return self.NT2VEC.get(code, self.NNTVEC)

    def _pr_encode_nt(self, a1_code, a2_code):
        """Encode allele pairs into matrix"""
        return self._pr_nt2vc(a1_code) + self._pr_nt2vc(a2_code)

    @staticmethod
    def _pr_try_open(file_path):
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

    @staticmethod
    def _pr_try_close(tables_obj):
        if tables_obj is None:
            return True

        try:
            tables_obj.close()
            return True
        except IOError as ioerr:
            print(ioerr, file=sys.stderr)
        except Exception as exp:
            print(exp, file=sys.stderr)

    def _pr_gen_gnm_region(self, gene_id, featuretype="transcript", which_mrna=None):
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

        # Using list to store each mRNA if necessary.
        self._CHROM = mrna.seqid
        self.mrna_pool[gene_id] = mrna
        self.exon_pool[gene_id] = exon_pool

    def _pr_gen_seq_mtrx(self, seq_itvl=None, shift=1e4):
        """Generate a chain of amb for a sequence.
        """
        seq_code_pool, snp_indx_pool, snp_code_pool, hap_code_pool, hap_phase_pool = self._pr_get_nodes(self._CHROM, ( "seq", "snp", "hap"))
        itvl_start, itvl_stop = int(seq_itvl.start - shift), seq_itvl.start

        sample_idx = self._pr_sample_id_to_idx(self._CHROM, self._SAMPLE_ID)
        _ntmtrx_pool = [None, None]
        seq_code = seq_code_pool[itvl_start: itvl_stop]
        snp_code = snp_code_pool[itvl_start: itvl_stop]
        ntmtrx = []
        for a1_code, a2_code in list(zip(seq_code, snp_code)):  # Scan the sequence one by one, could be too slow?
            if a2_code != -1 and hap_phase_pool[a2_code, sample_idx] == 1:  # Is phased
                snp_info = snp_indx_pool[a2_code]
                a1_idx, a2_idx = hap_code_pool[a2_code, sample_idx * 2: sample_idx * 2 + 2]
                a1_code, a2_code = snp_info[a1_idx + 2], snp_info[a2_idx + 2]
            else:
                a2_code = a1_code
            ntmtrx.append(self._pr_encode_nt(a1_code, a2_code))

        return np.array([ntmtrx])

    def _pr_gen_rcp(self, exon_pool):
        """Generate a ReadCountsPool.
        """
        ref_rc_pool, alt_rc_pool, snp_code_pool, snp_indx_pool, hap_code_pool, hap_phase_pool = self._pr_get_nodes(self._CHROM, ("rc", "snp", "hap"))

        sample_idx = self._pr_sample_id_to_idx(self._CHROM, self._SAMPLE_ID)
        rcp = ReadCountPool()
        for exon_itvl in exon_pool:
            exon_start, exon_stop = exon_itvl.start - 1, exon_itvl.stop - 1
            _id = (exon_itvl.seqid, exon_start, exon_stop)

            het_loci = snp_indx_pool[exon_start: exon_stop]
            ref_rc = ref_rc_pool[exon_start: exon_stop]
            alt_rc = alt_rc_pool[exon_start: exon_stop]

            has_het_loci = np.select(het_loci > -1, het_loci)  # FIXME: could be problematic.
            has_ref_rc = np.select(ref_rc > 0, ref_rc)
            has_alt_rc = np.select(alt_rc > 0, alt_rc)

            if has_het_loci and (has_ref_rc or has_alt_rc):
                nz_a1_counts, nz_a2_counts, nz_het_loci = [], [], []
                for rc, ac, loci_index in zip(ref_rc, alt_rc, het_loci):
                    if rc > 0 or ac > 0:
                        hap_code = hap_code_pool[loci_index, sample_idx * 2: sample_idx * 2 + 2]

                        if hap_code[0] == 1:
                            nz_a1_counts.append(rc)
                            nz_a2_counts.append(ac)
                        else:
                            nz_a1_counts.append(ac)
                            nz_a2_counts.append(rc)

                        nz_het_loci.append([snp_code_pool[loci_index], hap_code])

                rcp.add_a1_rc(_id, nz_a1_counts)
                rcp.add_a2_rc(_id, nz_a2_counts)
                rcp.add_loci_info(_id, nz_het_loci)

        return rcp

    def _pr_gen_ase_effect(self, exon_pool=None, mthd="bn", meta_exon=False):
        """Generate allele-specific expression effects from read counts.
        """
        rcp = self._pr_gen_rcp(exon_pool)
        if len(rcp):
            aefact = ASEeffectFactory(rcp)
            if mthd == "bb":
                return aefact.bbtest(), aefact.rc
            elif mthd == "bn":
                return aefact.bbtest(), aefact.rc

        return None

    def init(self, ant_db_name=None):
        """Initialize a factory by loading all needed files.
        """
        args = self.args
        self.hap_tab = self._pr_try_open(args.haplotypes)
        self.snp_idx = self._pr_try_open(args.snp_index)
        self.snp_tab = self._pr_try_open(args.snp_tab)
        self.seq_tab = self._pr_try_open(args.seq_tab)
        self.ref_tab = self._pr_try_open(args.ref_tab)
        self.alt_tab = self._pr_try_open(args.alt_tab)
        self._SAMPLE_ID = args.sample_id

        if ant_db_name is None:
            ant_db_name = os.path.splitext(args.genome_annot)[0] + ".db"

        if not os.path.exists(ant_db_name):
            gffutils.create_db(args.genome_annot, ant_db_name, disable_infer_transcripts=True, disable_infer_genes=True)

        self.ant_sql = gffutils.FeatureDB(ant_db_name)

        return self

    def gen_gnm_itvl(self, **kwargs):
        """Generate genomic intervals from Ensembl gene IDs.
        """
        self._parse_gene_ids()
        for _id in self._gene_ids:
            self._pr_gen_gnm_region(_id, **kwargs)

        return self

    def gen_seq_mtrx(self, shift_factor=1e4):
        """Generage regulatory sequence matrix.

        TODO:
            1. Add the parameters into logging file
        """
        shift = 2 * math.ceil(math.sqrt(shift_factor)) ** 2

        if self.mrna_pool:  # Not empty
            self.ntmtrx_pool = {gene_id: self._pr_gen_seq_mtrx(seq_itvl=mrna, shift=shift) for gene_id, mrna in self.mrna_pool.items()}
        else:
            print("It looks the self.mrna_pool is empty, did you use gen_gnm_region() yet?", file=sys.stderr)

        return self

    def gen_ase_effect(self, mthd="bn", meta_exon=False):
        if self.exon_pool:  # Not empty
            self.ase_pool = {
                gene_id: self._pr_gen_ase_effect(exon_pool=exon_pool, mthd=mthd, meta_exon=meta_exon)
                for gene_id, exon_pool in self.exon_pool.items()
            }
        else:
            print("It looks the self.exon_pool is empty, did you use gen_gnm_region() yet?", file=sys.stderr)

        return self

    def save_ase_report(self, opt_file=None):
        """Save the estimated ASE effects into disk.
        """
        if opt_file:
            opt_file = opt_file
        elif self.args.as_ase_report:
            opt_file = self.args.as_ase_report
        else:
            opt_file = self._SAMPLE_ID + "_ase_report.txt"

        header = "\t".join(["sample_id", "gene_id", "llh_ratio", "p_val", "ase_direction", "snp_ids", "snp_info", "snp_phase", "allele_counts"])
        with open(opt_file, "wt") as rfh:
            rfh.write(header + "\n")

            for gene_id, ase_effect in self.ase_pool.items():
                if ase_effect:
                    (llr, p_val, direction), (a1_rc, a2_rc, snp_info) = ase_effect
                else:
                    print("Skipping {} for no ASE effect information".format(gene_id))
                    continue

                snp_id_chain, snp_rc_chain, snp_pos_chain, snp_phase_chain = [""], ["A1|A2:"], ["CH,PS,REF,ALT:"], ["A1|A2:"]

                for _a1_rc, _a2_rc, ((snp_id, snp_pos, ref, alt), (a1_gntp, a2_gntp)) in zip(a1_rc, a2_rc, snp_info):
                    snp_id = snp_id.decode()
                    ref = ref.decode()
                    alt = alt.decode()

                    snp_pos = self._CHROM + "," + str(snp_pos) + "," + ref + ">" + alt
                    snp_phase = str(a1_gntp) + "|" + str(a2_gntp)
                    snp_rc = str(_a1_rc) + "|" + str(_a2_rc)

                    if snp_id == ".":
                        snp_id = snp_pos

                    snp_id_chain.append(snp_id)
                    snp_rc_chain.append(snp_rc)
                    snp_pos_chain.append(snp_pos)
                    snp_phase_chain.append(snp_phase)

                info_list_1 = [self._SAMPLE_ID, gene_id, str(llr), str(p_val), str(direction)]
                info_list_2 = [x[0] + ";".join(x[1:]) for x in [snp_id_chain, snp_pos_chain, snp_phase_chain, snp_rc_chain]]

                est_result = "\t".join(info_list_1 + info_list_2)
                rfh.write(est_result + "\n")

        return self

    def save_train_set(self, opt_file=None):
        if opt_file:
            opt_file = opt_file
        elif self.args.as_train_set:
            opt_file = self.args.as_train_set
        else:
            opt_file = self._SAMPLE_ID + "_matrix_and_ase.npz"

        output_dataset = {}
        for (gene_id, ntmtrx), (_, ase) in zip(self.ntmtrx_pool.items(), self.ase_pool.items()):
            if ase:
                _ase_effect = np.array([ntmtrx, list(ase[0][1:])])
            else:
                _ase_effect = np.array([ntmtrx, [1, 0]])

            output_dataset[gene_id] = _ase_effect

        np.savez(opt_file, **output_dataset, allow_pickle=True)

        return self

    def shutdown(self):
        self._pr_try_close(self.hap_tab)
        self._pr_try_close(self.snp_idx)
        self._pr_try_close(self.snp_tab)
        self._pr_try_close(self.ref_tab)
        self._pr_try_close(self.alt_tab)
        self._pr_try_close(self.seq_tab)

if __name__ == "__main__":
    print("[W]: This module should not be executed directly.", file=sys.stderr)
