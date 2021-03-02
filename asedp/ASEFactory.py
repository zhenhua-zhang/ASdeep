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

logging.basicConfig(format='{levelname: ^8}| {asctime} | {name} | {message}', style='{',
                    datefmt='%Y-%m-%d, %H:%M:%S', level=logging.INFO)

from collections import UserDict

import vcf
import pyfaidx
import gffutils

import pandas as pd
import numpy as np
from scipy.optimize import minimize
from scipy.stats import betabinom, binom, chi2

from zutils import cmp
from zutils import M2B

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
        for _, _value in self.data.items():
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

        # _est_t = opt_results["x"]
        estll = opt_results["fun"]

        hyp_t = 0.5
        hypll = self._pr_bnllh(hyp_t, k_vec, n_vec)

        llr = - 2 * (estll - hypll)
        p_val = chi2.sf(llr, 1)

        k_sum, n_sum = int(sum(k_vec)), int(sum(n_vec))
        return llr, p_val, cmp(2 * k_sum, n_sum) # likelihood ratio, p-value, direction

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

        # _est_a, _est_b = opt_results["x"]
        estll = opt_results["fun"]

        hyp_a = hyp_b = opt_results["x"].mean()
        hypll = self._pr_bbllh([hyp_a, hyp_b], k_vec, n_vec)

        llr = - 2 * (estll - hypll) # The likelihood value is timed by -1
        p_val = chi2.sf(llr, 1)

        k_sum, n_sum = int(sum(k_vec)), int(sum(n_vec))
        return llr, p_val, cmp(2 * k_sum, n_sum) # likelihood ratio, p-value, direction

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
    _sample_id = None  # Static
    _gene_ids = None   # Static

    variant_pool = None
    sequence_pool = None
    readcount_pool = None
    gene_features = None

    def __init__(self, args):
        super(ASEFactory, self).__init__()
        self.args = args

        self._sample_id = args.sample_id
        self.snp_pool = {}
        self.mrna_pool = {}
        self.exon_pool = {}
        self.ase_pool = None
        self.seq_pool = None

        if self.args.variants_file:
            self.variant_pool = vcf.Reader(open(self.args.variants_file, 'rb'))
        else:
            raise FileNotFoundError("No variants file was given or it's not accessible!")

        if self.args.genome_seq_file:
            self.sequence_pool = pyfaidx.Fasta(self.args.genome_seq_file)
        else:
            raise FileNotFoundError("No genome sequence file was given or it's not accessible!")

        if self.args.ase_readcount_file:
            self.ase_readcounts = pd.read_csv(self.args.ase_readcount_file, sep=' ')
        else:
            raise FileNotFoundError("No ASE read counts file was given or it's not accessible!")

        if self.args.gene_feature_file:
            gene_feature_file = self.args.gene_feature_file
            ant_db_name = os.path.splitext(gene_feature_file)[0] + ".db"
            if not os.path.exists(ant_db_name):
                gffutils.create_db(args.gene_feature_file, ant_db_name,
                                   disable_infer_transcripts=True, disable_infer_genes=True)
            self.gene_features = gffutils.FeatureDB(ant_db_name)
        else:
            raise FileNotFoundError("No gene feature file was given or it's not accessible!")

        self._parse_gene_ids()

    def _parse_gene_ids(self):
        self._gene_ids = self.args.gene_ids + self._parse_gene_id_file()

    def _parse_gene_id_file(self):
        if self.args.gene_id_file:
            with open(self.args.gene_id_file) as gifh:
                ids = [x.strip("\n") for x in gifh]
            return ids
        return []
 
    def _pr_get_genome_itvl(self, gene_id, featuretype="transcript", which_mrna=None):
        # Get genomic interval for the given gene_id
        mrna_pool = self.gene_features.children(gene_id, featuretype=featuretype)

        if which_mrna is None:
            mrna = next(mrna_pool, False)
        else:
            raise NotImplementedError("Not implemented yet")

        if mrna:
            parent_id = mrna.id
        else:
            parent_id = gene_id

        exon_pool = self.gene_features.children(parent_id, featuretype="exon")

        # Using list to store each mRNA if necessary.
        self.mrna_pool[gene_id] = mrna
        self.exon_pool[gene_id] = exon_pool

    def _pr_get_variants(self, chrom, start, end) -> vcf.Reader:
        # Get variants in the genomic interval
        # 0-based
        return self.variant_pool.fetch(chrom, start, end)

    def _pr_get_sequence(self, chrom, start, end) -> pyfaidx.Fasta:
        # Get sequence in the genomic interval
        # 1-based
        return self.sequence_pool.get_seq(chrom, start, end)

    def _pr_get_readcounts(self, chrom, start, end) -> pd.DataFrame:
        # Get the read counts in the genomic interval
        # TODO: check the WASP manual to determine x-based
        # pdb.set_trace()
        return self.ase_readcounts.query('chrom==@chrom & @start <= pos & pos <= @end')

    @staticmethod
    def _pr_encode_bnt(gntp):
        return M2B.get(gntp, "N")

    def _pr_gen_ase_seq(self, itvl: gffutils.Feature, shift=5e2):
        #pdb.set_trace()
        chrom, start, end, strand = itvl.chrom, itvl.start, itvl.end, itvl.strand
        if strand == '+':
            start, end = start - shift, start
        elif strand == '-':
            start, end = end, end + shift
        else:
            logging.warning('No strand found, skip interval: {}:{}-{}'.format(chrom, start, end))
            return None

        sequence = self._pr_get_sequence(chrom, start, end)
        variants = self._pr_get_variants(chrom, start, end)

        amb_sequence = sequence.seq
        for vcf_record in variants:
            sample_genotype = vcf_record.genotype(self._sample_id)

            if sample_genotype.is_het and sample_genotype.phased:
                genotype_bases = sample_genotype.gt_bases.replace('|', '')
                amb_base = self._pr_encode_bnt(genotype_bases)

                insert_pos = vcf_record.start - start
                amb_sequence = sequence[:insert_pos].seq + amb_base + sequence[insert_pos+1:].seq
            else:
                logging.warning('Variant not phased, skip: {}:{}'.format(chrom, vcf_record.end))

        return amb_sequence

    # TODO: rewrite this function
    def _pr_gen_rcp(self, exon_pool):
        # Generate a ReadCountsPool.
        # TODO: check x-based
        rcp = ReadCountPool()
        for exon_itvl in exon_pool:
            exon_chrom, exon_start, exon_stop = exon_itvl.seqid, exon_itvl.start - 1, exon_itvl.stop
            _id = (exon_chrom, exon_start, exon_stop)

            exon_vars = self._pr_get_readcounts(exon_chrom, exon_start, exon_stop)
            if exon_vars.empty:
                continue

            het_exon_vars = exon_vars.loc[exon_vars.apply(lambda x: x['gt'] in ['0|1', '1|0'], axis=1), :]
            if het_exon_vars.empty:
                continue

            # pdb.set_trace()
            nz_a1_counts, nz_a2_counts, nz_het_loci = [], [], []
            for _, (chrom, pos, ref, alt, gt, ref_rc, alt_rc, _) in het_exon_vars.iterrows():
                if gt == '0|1':
                    nz_a1_counts.append(ref_rc)
                    nz_a2_counts.append(alt_rc)
                else:
                    nz_a1_counts.append(alt_rc)
                    nz_a2_counts.append(ref_rc)

                nz_het_loci.append([[chrom, pos, ref, alt], [gt[0], gt[-1]]])

            rcp.add_a1_rc(_id, nz_a1_counts)
            rcp.add_a2_rc(_id, nz_a2_counts)
            rcp.add_loci_info(_id, nz_het_loci)

        return rcp

    def _pr_gen_ase(self, exon_pool=None, mthd="bn"):
        """Generate allele-specific expression effects from read counts.
        """
        rcp = self._pr_gen_rcp(exon_pool)
        # pdb.set_trace()
        if len(rcp) != 0:
            aefact = ASEeffectFactory(rcp)
            if mthd == "bb":
                return aefact.bbtest(), aefact.read_counts
            elif mthd == "bn":
                return aefact.bntest(), aefact.read_counts

        return None


    def gen_gnm_itvl(self, **kwargs):
        """Generate genomic intervals from Ensembl gene IDs.
        """
        for _id in self._gene_ids:
            self._pr_get_genome_itvl(_id, **kwargs)

        return self

    def gen_seq(self, shift_factor=1e3):
        """Generage regulatory sequence matrix.
        """
        shift = 2 * math.ceil(math.sqrt(shift_factor)) ** 2

        if self.mrna_pool:  # Not empty
            self.seq_pool = {
                gene_id: self._pr_gen_ase_seq(itvl=mrna, shift=shift)
                for gene_id, mrna in self.mrna_pool.items()}
        else:
            logging.warning("The mrna_pool is empty, use gen_gnm_itvl() first.")

        return self

    def gen_ase(self, mthd="bn"):
        """Generate ASE effects.
        """
        if self.exon_pool:  # Not empty
            self.ase_pool = {
                gene_id: self._pr_gen_ase(exon_pool, mthd)
                for gene_id, exon_pool in self.exon_pool.items()}
        else:
            logging.warn("The exon_pool is empty, use gen_gnm_itvl() first.")

        return self

    def save_ase_report(self, opt_file=None):
        """Save the estimated ASE effects into disk.
        """
        if opt_file is None:
            if self.args.as_ase_report:
                opt_file = self.args.as_ase_report
            else:
                opt_file = self._sample_id + ".ase_report.txt"

        header = "\t".join(["sample_id", "gene_id", "llh_ratio", "p_val", "ase_direction",
                            "snp_info", "snp_phase", "allele_counts"])
        with open(opt_file, "wt") as rfh:
            rfh.write(header + "\n")

            for gene_id, ase_effect in self.ase_pool.items():
                if ase_effect:
                    (llr, p_val, direction), (a1_rc, a2_rc, snp_info) = ase_effect
                    snp_rc_chain, snp_pos_chain, snp_phase_chain = ["A1|A2:"], ["CH,PS,REF>ALT:"], ["A1|A2:"]
                    # pdb.set_trace()
                    for _a1_rc, _a2_rc, ((chrom, pos, ref, alt), (a1_gntp, a2_gntp)) in zip(a1_rc, a2_rc, snp_info):
                        snp_pos = str(chrom) + "," + str(pos) + "," + ref + ">" + alt
                        snp_phase = str(a1_gntp) + "|" + str(a2_gntp)
                        snp_rc = str(_a1_rc) + "|" + str(_a2_rc)

                        snp_rc_chain.append(snp_rc)
                        snp_pos_chain.append(snp_pos)
                        snp_phase_chain.append(snp_phase)

                    info_list_1 = [self._sample_id, gene_id, str(llr), str(p_val), str(direction)]
                    info_list_2 = [x[0] + ";".join(x[1:]) for x in [snp_pos_chain, snp_phase_chain, snp_rc_chain]]
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
        for _gene_id in self.seq_pool.keys():
            _ntseq = self.seq_pool[_gene_id]
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


if __name__ == "__main__":
    print("[W]: This module should not be executed directly.", file=sys.stderr)
