"""ASE factory version 0.2.0

TODO:
    1. In output FASTA file, only gene id is used for index. However, there
    could be duplicated records for one gene id. Therefore, mRNA id for the gene
    id should be used instead of gene id only.
"""

import textwrap
import logging
import gzip
import math
import pdb
import sys
import os

from collections import UserDict

import numpy as np

import gffutils
import tables
from scipy.optimize import minimize
from scipy.stats import betabinom, binom, chi2

from zutils import cmp, logger
from zutils import M2B

NNTVEC = [0, 0, 0, 0]  # Base not in ACGTN
NT2VEC = {"A": [1, 0, 0, 0], "C": [0, 1, 0, 0], "G": [0, 0, 1, 0], "T": [0, 0, 0, 1]}

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
        return default

    def __len__(self):
        return len(self.data)

    def _pr_add(self, itvl_id: str, content: list, add_type="a1"):
        """Add counts.
        """
        type_pool = ["a1", "a2", "info"]
        if itvl_id not in self.data:
            self.data[itvl_id] = [[], [], []]

        self.data[itvl_id][type_pool.index(add_type)] = content

    def add_a1_rc(self, itvl_id, counts):
        self._pr_add(itvl_id, counts, add_type="a1")

    def add_a2_rc(self, itvl_id, counts):
        self._pr_add(itvl_id, counts, add_type="a2")

    def add_loci_info(self, itvl_id, info):
        self._pr_add(itvl_id, info, add_type="info")

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

    def __init__(self, read_counts):
        self.read_counts = self._pr_parse_rc(read_counts)
        self.optim_results = {}

    @staticmethod
    def _pr_parse_rc(read_counts):
        if isinstance(read_counts, (tuple, list)) and len(read_counts) >= 2:
            return read_counts, []
        elif isinstance(read_counts, ReadCountPool):
            return read_counts.unpack()
        else:
            raise TypeError("read_counts parameter should be tuple, list, or ReadCountPool")

    @staticmethod
    def _pr_bnllh(param, k_vec, n_vec):
        """The likelihood function of Binomial distribution.
        """
        lh_vec = [binom.logpmf(k, n, param) for k, n in zip(k_vec, n_vec)]
        return -sum(lh_vec)

    def bntest(self):
        """Do Binomial test on given data.
        """
        x_zero = 0.5
        k_vec = self._pr_get_k_vec()
        n_vec = self._pr_get_n_vec()

        opt_results = minimize(
            self._pr_bnllh, x_zero, args=(k_vec, n_vec),
            method="nelder-mead", options={"maxfev": 1e3, "ftol": 1e-8}
        )
        self.optim_results["binomial"] = opt_results

        _est_t = opt_results["x"]
        estll = opt_results["fun"]

        hyp_t = 0.5
        hypll = self._pr_bnllh(hyp_t, k_vec, n_vec)

        llr = - 2 * (estll - hypll)
        p_val = chi2.sf(llr, 1)

        k_sum, n_sum = int(sum(k_vec)), int(sum(n_vec))
        return llr, p_val, cmp(2 * k_sum, n_sum)  # likelihood ratio, p-value, direction

    @staticmethod
    def _pr_bbllh(params, k_vec, n_vec):
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
        x_zero = [0.5, 0.5]
        k_vec = self._pr_get_k_vec()
        n_vec = self._pr_get_n_vec()

        opt_results = minimize(
            self._pr_bbllh, x_zero, args=(k_vec, n_vec),
            method="nelder-mead", options={"maxfev": 1e3, "ftol": 1e-8}
        )
        self.optim_results["beta_binomial"] = opt_results

        _est_a, _est_b = opt_results["x"]
        estll = opt_results["fun"]

        hyp_a = hyp_b = opt_results["x"].mean()
        hypll = self._pr_bbllh([hyp_a, hyp_b], k_vec, n_vec)

        llr = - 2 * (estll - hypll)  # The likelihood value is timed by -1
        p_val = chi2.sf(llr, 1)

        k_sum, n_sum = int(sum(k_vec)), int(sum(n_vec))
        return llr, p_val, cmp(2 * k_sum, n_sum)  # likelihood ratio, p-value, direction

    def _pr_get_k_vec(self):
        return self.read_counts[0]

    def _pr_get_n_vec(self):
        return list(map(sum, zip(self.read_counts[0], self.read_counts[1])))

    def chi_square_test(self):
        """Using a Chi-square test on given data.
        """
        return NotImplemented


class ASEFactory:
    """A class to produce ASE effects and matrix of regulatory sequence.
    """
    _sample_id_to_idx = None  # Has to be None for there's a if-else branch
    _chrom = None  # Static
    _sample_id = None  # Static
    _gene_ids = None
    hap_tab, snp_idx, snp_tab, seq_tab, ref_tab, alt_tab, ant_sql = [None] * 7

    def __init__(self, args):
        super(ASEFactory, self).__init__()
        self.args = args
        self.mrna_pool = {}
        self.exon_pool = {}
        self.ase_pool = None
        self.ntseq_pool = None

    def _parse_gene_ids(self):
        self._gene_ids = self.args.gene_ids + self._parse_gene_id_file()

    def _parse_gene_id_file(self):
        if self.args.gene_id_file:
            with open(self.args.gene_id_file) as gifh:
                ids = [x.strip("\n") for x in gifh]
            return ids
        return []

    def _pr_sample_id_to_idx(self, chrom, sample_id):
        if self._sample_id_to_idx is None:
            hap_samples = self.hap_tab.get_node("/samples_{}".format(chrom))
            self._sample_id_to_idx = {
                _sample[0].decode(): _idx
                for _idx, _sample in enumerate(hap_samples)
            }

        return self._sample_id_to_idx.get(sample_id)

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
                logger.warn("Unsupported type of node: " + node_name)

        return tuple(node_pool)

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
        self._chrom = mrna.seqid
        self.mrna_pool[gene_id] = mrna
        self.exon_pool[gene_id] = exon_pool

    def _pr_encode_bnt(self, a1_code, a2_code):
        return M2B.get((a1_code, a2_code), "N")

    def _pr_gen_seq(self, seq_itvl=None, shift=5e2):
        # pdb.set_trace()
        seq_code_pool, snp_indx_pool, snp_code_pool, hap_code_pool, hap_phase_pool = self._pr_get_nodes(self._chrom, ( "seq", "snp", "hap"))
        itvl_start, itvl_stop = int(seq_itvl.start - shift), seq_itvl.start

        sample_idx = self._pr_sample_id_to_idx(self._chrom, self._sample_id)
        seq_code = seq_code_pool[itvl_start: itvl_stop]
        snp_code = snp_code_pool[itvl_start: itvl_stop]

        ntstr = ""
        # Scan the sequence one by one, could be too slow?
        for a1_code, a2_code in list(zip(seq_code, snp_code)):
            if a2_code != -1 and hap_phase_pool[a2_code, sample_idx] == 1:
                snp_info = snp_indx_pool[a2_code]
                a1_idx, a2_idx = hap_code_pool[a2_code, sample_idx * 2: sample_idx * 2 + 2]
                if a1_idx == -1 or a2_idx == -1:
                    a2_code = a1_code
                else:
                    a1_code, a2_code = snp_info[a1_idx + 2][0], snp_info[a2_idx + 2][0]
            else:
                a2_code = a1_code

            ntstr += self._pr_encode_bnt(a1_code, a2_code)

        return ntstr

    def _pr_gen_rcp(self, exon_pool):
        """Generate a ReadCountsPool.
        """
        # pdb.set_trace()
        ref_rc_pool, alt_rc_pool, snp_code_pool, snp_indx_pool, hap_code_pool, hap_phase_pool = self._pr_get_nodes(self._chrom, ("rc", "snp", "hap"))

        sample_idx = self._pr_sample_id_to_idx(self._chrom, self._sample_id)
        rcp = ReadCountPool()
        for exon_itvl in exon_pool:
            exon_start, exon_stop = exon_itvl.start - 1, exon_itvl.stop - 1
            _id = (exon_itvl.seqid, exon_start, exon_stop)

            het_loci = snp_indx_pool[exon_start: exon_stop]
            ref_rc = ref_rc_pool[exon_start: exon_stop]
            alt_rc = alt_rc_pool[exon_start: exon_stop]

            has_het_loci = np.select(het_loci > -1, het_loci)
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

    def _pr_gen_ase(self, exon_pool=None, mthd="bn", meta_exon=False):
        """Generate allele-specific expression effects from read counts.
        """
        rcp = self._pr_gen_rcp(exon_pool)
        # pdb.set_trace()
        if len(rcp) != 0:
            aefact = ASEeffectFactory(rcp)
            if mthd == "bb":
                return aefact.bbtest(), aefact.read_counts
            elif mthd == "bn":
                return aefact.bbtest(), aefact.read_counts

        return None

    def init(self, ant_db_name=None):
        """Initialize a factory by loading all needed files.

        Arguments:
            ant_db_name (None, str): The path to gene feature format annotations
            database converted by package `gffutils`.
        """
        args = self.args
        self.hap_tab = self._pr_try_open(args.haplotypes)
        self.snp_idx = self._pr_try_open(args.snp_index)
        self.snp_tab = self._pr_try_open(args.snp_tab)
        self.seq_tab = self._pr_try_open(args.seq_tab)
        self.ref_tab = self._pr_try_open(args.ref_tab)
        self.alt_tab = self._pr_try_open(args.alt_tab)
        self._sample_id = args.sample_id

        if ant_db_name is None:
            ant_db_name = os.path.splitext(args.genome_annot)[0] + ".db"

        if not os.path.exists(ant_db_name):
            gffutils.create_db(args.genome_annot, ant_db_name,
                               disable_infer_transcripts=True,
                               disable_infer_genes=True)

        self.ant_sql = gffutils.FeatureDB(ant_db_name)
        self._parse_gene_ids()

        return self

    def gen_gnm_itvl(self, **kwargs):
        """Generate genomic intervals from Ensembl gene IDs.
        """
        for _id in self._gene_ids:
            self._pr_gen_gnm_region(_id, **kwargs)

        return self

    def gen_seq(self, shift_factor=1e3):
        """Generage regulatory sequence matrix.
        """
        shift = 2 * math.ceil(math.sqrt(shift_factor)) ** 2

        if self.mrna_pool:  # Not empty
            self.ntseq_pool = {
                gene_id: self._pr_gen_seq(seq_itvl=mrna, shift=shift)
                for gene_id, mrna in self.mrna_pool.items()}
        else:
            logger.warning("The mrna_pool is empty, use gen_gnm_region() first.")

        return self

    def gen_ase(self, mthd="bn", meta_exon=False):
        """Generate ASE effects.
        """
        if self.exon_pool:  # Not empty
            self.ase_pool = {
                gene_id: self._pr_gen_ase(exon_pool, mthd, meta_exon)
                for gene_id, exon_pool in self.exon_pool.items()}
        else:
            logger.warn("The exon_pool is empty, use gen_gnm_region() first.")

        return self

    def save_ase_report(self, opt_file=None):
        """Save the estimated ASE effects into disk.
        """
        if opt_file is None:
            if self.args.as_ase_report:
                opt_file = self.args.as_ase_report
            else:
                opt_file = self._sample_id + ".ase_report.txt"

        header = "\t".join(["sample_id", "gene_id", "llh_ratio", "p_val", "ase_direction", "snp_ids", "snp_info", "snp_phase", "allele_counts"])
        with open(opt_file, "wt") as rfh:
            rfh.write(header + "\n")

            for gene_id, ase_effect in self.ase_pool.items():
                if ase_effect:
                    (llr, p_val, direction), (a1_rc, a2_rc, snp_info) = ase_effect
                    snp_id_chain, snp_rc_chain, snp_pos_chain, snp_phase_chain = [""], ["A1|A2:"], ["CH,PS,REF,ALT:"], ["A1|A2:"]
                    for _a1_rc, _a2_rc, ((snp_id, snp_pos, ref, alt), (a1_gntp, a2_gntp)) in zip(a1_rc, a2_rc, snp_info):
                        snp_id, ref, alt = snp_id.decode(), ref.decode(), alt.decode()

                        snp_pos = self._chrom + "," + str(snp_pos) + "," + ref + ">" + alt
                        snp_phase = str(a1_gntp) + "|" + str(a2_gntp)
                        snp_rc = str(_a1_rc) + "|" + str(_a2_rc)

                        if snp_id == ".":
                            snp_id = snp_pos

                        snp_id_chain.append(snp_id)
                        snp_rc_chain.append(snp_rc)
                        snp_pos_chain.append(snp_pos)
                        snp_phase_chain.append(snp_phase)

                    info_list_1 = [self._sample_id, gene_id, str(llr), str(p_val), str(direction)]
                    info_list_2 = [x[0] + ";".join(x[1:]) for x in [snp_id_chain, snp_pos_chain, snp_phase_chain, snp_rc_chain]]
                else:
                    info_list_1 = [self._sample_id, gene_id, "NA", "NA", "NA"]
                    info_list_2 = ["NA"] * 4

                est_result = "\t".join(info_list_1 + info_list_2)
                rfh.write(est_result + "\n")

        return self

    def save_train_set(self, opt_file=None, save_fmt="fa.gz"):
        """Save seqeunce and ASE effects into disk.

        Args:
            opt_file (str, optional, None): The output file.
            save_fmt (str, optional, fa.gz): output format. ["fa", "fa.gz"]

        Returns:
            self (:obj:`ASEFactory`): The object itself.

        Raises:
            ValueError: If save_fmt argument is not one of [fa, fa.gz]
        """
        if opt_file is None:
            if self.args.as_train_set:
                opt_file = self.args.as_train_set
            else:
                opt_file = self._sample_id + ".ntsq_and_ase." + save_fmt

        _opt_str = ""
        _sample_id = self._sample_id
        for _gene_id in self.ntseq_pool.keys():
            _ntseq = self.ntseq_pool[_gene_id]
            _ase = self.ase_pool[_gene_id]

            _ase_effect = "{}|{}".format(*list(_ase[0][1:])) if _ase else "1|0"
            _nt_seq_fold = "\n".join(textwrap.wrap(_ntseq))

            _opt_str += ">{}|{}|{}\n{}\n".format(_sample_id, _gene_id, _ase_effect, _nt_seq_fold)

        if save_fmt == "fa":
            _open = open
        elif save_fmt == "fa.gz":
            _open = gzip.open
            _opt_str = _opt_str.encode()
        else:
            raise TypeError("Unknown type of output format: {}".format(save_fmt))

        with _open(opt_file, "w") as _opfh:
            _opfh.write(_opt_str)

        return self

    def shutdown(self):
        """Close all open HDF5 files.
        """
        self._pr_try_close(self.hap_tab)
        self._pr_try_close(self.snp_idx)
        self._pr_try_close(self.snp_tab)
        self._pr_try_close(self.ref_tab)
        self._pr_try_close(self.alt_tab)
        self._pr_try_close(self.seq_tab)


if __name__ == "__main__":
    print("[W]: This module should not be executed directly.", file=sys.stderr)
