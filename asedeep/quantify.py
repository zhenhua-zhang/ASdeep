"""ASE factory version 0.2.0

TODO:
    1. In output FASTA file, only gene id is used for index. However, there
    could be duplicated records for one gene id. Therefore, mRNA id for the gene
    id should be used instead of gene id only.
"""

# Built-in packages
import os
import math
import gzip
import logging
import textwrap

# Third-party packages
import vcf
import pandas
import pyfaidx
import gffutils

from scipy.optimize import minimize
from scipy.stats import betabinom, binom, chi2

# Modules from this package
from .zutils import cmp, M2B

class ReadCountPool:
    """A class for Read counts."""
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

    # Add counts
    def _pr_add(self, itvl_id: str, content: list, add_type="a1"):
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
        """Unpack the object."""
        a1_rc_pool, a2_rc_pool, info_pool = [], [], []
        for _, _value in self.data.items():
            a1_rc_pool.extend(_value[0])
            a2_rc_pool.extend(_value[1])
            info_pool.extend(_value[2])
        return a1_rc_pool, a2_rc_pool, info_pool


class ASEeffectFactory:
    """Estimate the allele specific expression effects using likelihood ratio test under
    Binomial / Beta-Binomial distribution.

    TODO:
        1. Not able to deal with meta-exons.
        2. Which transcript to be considered.
    """

    def __init__(self, read_counts, min_ac=3, min_gc=5):
        self.read_counts = self._pr_parse_rc(read_counts)
        self.optim_results = {}
        self.min_gc = min_gc
        self.min_ac = min_ac

    @staticmethod
    def _pr_parse_rc(read_counts):
        if isinstance(read_counts, (tuple, list)) and len(read_counts) == 2:
            return read_counts, []
        elif isinstance(read_counts, ReadCountPool):
            return read_counts.unpack()
        else:
            raise TypeError("read_counts parameter should be tuple, list, or ReadCountPool")

    # The likelihood function of Binomial distribution
    @staticmethod
    def _pr_bnllh(param, k_vec, n_vec):
        lh_vec = [binom.logpmf(k, n, param) for k, n in zip(k_vec, n_vec)]
        return -sum(lh_vec)

    def bntest(self):
        """Do Binomial test on given data."""
        x_zero = 0.5
        k_vec = self._pr_get_k_vec()
        n_vec = self._pr_get_n_vec()

        a1c = int(sum(k_vec))
        aac = int(sum(n_vec))
        a2c = aac - a1c

        do_ase = (a1c >= self.min_ac or a2c >= self.min_ac) and aac >= self.min_gc
        if do_ase:
            opt_results = minimize(self._pr_bnllh, x_zero, args=(k_vec, n_vec),
                                   method="nelder-mead", options={"maxfev": 1e3, "ftol": 1e-8})
            self.optim_results["binomial"] = opt_results

            # _est_t = opt_results["x"]
            estll = opt_results["fun"]

            hyp_t = 0.5
            hypll = self._pr_bnllh(hyp_t, k_vec, n_vec)

            llr = - 2 * (estll - hypll)
            p_val = chi2.sf(llr, 1)

            k_sum, n_sum = sum(k_vec), sum(n_vec)
            return llr, p_val, cmp(2 * k_sum, n_sum) # likelihood ratio, p-value, direction

        return 0, 1, 0

    # The likelihood function of Beta-Binomial distribution
    @staticmethod
    def _pr_bbllh(params, k_vec, n_vec):
        if len(params) != 2:
            raise TypeError("The params should be a list with two elements.")

        bb_a, bb_b = params
        lh_vec = [betabinom.logpmf(k, n, bb_a, bb_b) for k, n in zip(k_vec, n_vec)]
        return -sum(lh_vec)

    def bbtest(self):
        """Do Beta-binomial test on given data."""
        x_zero = [0.5, 0.5]
        k_vec = self._pr_get_k_vec()
        n_vec = self._pr_get_n_vec()

        a1c = int(sum(k_vec))
        aac = int(sum(n_vec))
        a2c = aac - a1c

        do_ase = (a1c >= self.min_ac or a2c >= self.min_ac) and aac >= self.min_gc
        if do_ase:
            opt_results = minimize(self._pr_bbllh, x_zero, args=(k_vec, n_vec),
                                   method="nelder-mead", options={"maxfev": 1e3, "ftol": 1e-8})
            self.optim_results["beta_binomial"] = opt_results

            # _est_a, _est_b = opt_results["x"]
            estll = opt_results["fun"]

            hyp_a = hyp_b = opt_results["x"].mean()
            hypll = self._pr_bbllh([hyp_a, hyp_b], k_vec, n_vec)

            llr = - 2 * (estll - hypll)
            p_val = chi2.sf(llr, 1)

            k_sum, n_sum = sum(k_vec), sum(n_vec)
            return llr, p_val, cmp(2 * k_sum, n_sum) # likelihood ratio, p-value, direction

        return 0, 1, 0

    def _pr_get_k_vec(self):
        return self.read_counts[0]

    def _pr_get_n_vec(self):
        return list(map(sum, zip(self.read_counts[0], self.read_counts[1])))

    def chi_square_test(self):
        """Using a Chi-square test on given data."""
        self.__str__
        return NotImplemented


class Quantifier:
    """A class to produce ASE effects and matrix of regulatory sequence."""
    _sample_id = None  # Static
    _gene_ids = None   # Static

    def __init__(self, args):
        super(Quantifier, self).__init__()
        self.args = args

        self._sample_id = args.sample_id
        self.snp_pool = {}
        self.mrna_pool = {}
        self.exon_pool = {}
        self.ase_pool = None
        self.seq_pool = None

        if self.args.variants_file:
            self.variant_pool = vcf.Reader(open(self.args.variants_file, "rb"))
        else:
            raise FileNotFoundError("No variants file was given or it's not accessible!")

        if self.args.genome_seq_file:
            self.sequence_pool = pyfaidx.Fasta(self.args.genome_seq_file)
        else:
            raise FileNotFoundError("No genome sequence file was given or it's not accessible!")

        if self.args.ase_readcount_file:
            self.ase_readcounts = pandas.read_csv(self.args.ase_readcount_file, sep=None,
                                                  engine="python", converters={"contig": str})
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

        self._gene_ids = self.args.gene_ids + self._parse_gene_id_file()

    # Gene id files
    def _parse_gene_id_file(self):
        if self.args.gene_id_file:
            with open(self.args.gene_id_file) as gifh:
                ids = [x.strip("\n") for x in gifh]
            return ids
        return []
 
    # Generate genomic intervals
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

    # Get variants in the genomic interval. 0-based
    def _pr_get_variants(self, chrom, start, end) -> vcf.Reader:
        return self.variant_pool.fetch(chrom, start-1, end)

    # Get sequence in the genomic interval. 1-based
    def _pr_get_sequence(self, chrom, start, end):
        return self.sequence_pool.get_seq(chrom, start, end)

    # Get the read counts in the genomic interval
    def _pr_get_readcounts(self, chrom, start, end) -> pandas.DataFrame:
        query_str = "contig == '{}' & {} <= position & position <= {}".format(chrom, start, end)
        return (self.ase_readcounts.query(query_str))

    # Fetch ambiguous
    @staticmethod
    def _pr_encode_bnt(gntp):
        return M2B.get(gntp, "N")

    # Generate allelic sequence using ambiguous base
    def _pr_gen_ase_seq(self, itvl: gffutils.Feature, shift=5e2):
        chrom, start, end, strand = itvl.chrom, itvl.start, itvl.end, itvl.strand
        if strand == "+":
            start, end = max(1, start - shift), start
        else:
            start, end = end, end + shift

        sequence = self._pr_get_sequence(chrom, start, end)
        variants = self._pr_get_variants(chrom, start, end)

        amb_sequence = sequence.seq
        for vcf_record in variants:
            sample_genotype = vcf_record.genotype(self._sample_id)

            if sample_genotype.is_het and sample_genotype.phased:
                genotype_bases = sample_genotype.gt_bases.replace("|", "")
                amb_base = self._pr_encode_bnt(genotype_bases)

                insert_pos = vcf_record.start - start
                amb_sequence = sequence[:insert_pos].seq + amb_base + sequence[insert_pos+1:].seq

        return amb_sequence

    # Generate a ReadCountsPool.
    def _pr_gen_rcp(self, exon_pool):
        def _get_var_gt(key, vmap):
            call = vmap.get((str(key.loc["contig"]), key.loc["position"]))
            if call:
                return call.data.GT

            return call

        rcp = ReadCountPool()
        for exon_itvl in exon_pool:
            exon_chrom, exon_start, exon_end = exon_itvl.seqid, exon_itvl.start, exon_itvl.stop

            exon_rc = self._pr_get_readcounts(exon_chrom, exon_start, exon_end)
            if exon_rc.empty:
                continue

            # Fetch variant information
            vmap = {(var.CHROM, var.POS): var.genotype(self._sample_id)
                    for var in self._pr_get_variants(exon_chrom, exon_start, exon_end)}
            # Add variant information to the allelic read count DataFrame
            exon_rc.insert(7, "genotype", exon_rc.apply(_get_var_gt, 1, vmap=vmap))
            # Only working on heterzygous loci
            het_exon_vars = exon_rc.loc[exon_rc["genotype"].apply(lambda x: x in ["0|1", "1|0"]), :]

            if het_exon_vars.empty:
                continue

            nz_a1_counts, nz_a2_counts, nz_het_loci = [], [], []
            for _, (chrom, pos, snpid, ref, alt, refrc, altrc, gt, *_) in het_exon_vars.iterrows():
                if gt == "0|1":
                    nz_a1_counts.append(refrc)
                    nz_a2_counts.append(altrc)
                else:
                    nz_a1_counts.append(altrc)
                    nz_a2_counts.append(refrc)

                if snpid == ".":
                    snpid = str(chrom) + ":" + str(pos)

                nz_het_loci.append([[chrom, pos, snpid, ref, alt], [gt[0], gt[-1]]])

            _id = (exon_chrom, exon_start, exon_end)
            rcp.add_a1_rc(_id, nz_a1_counts)
            rcp.add_a2_rc(_id, nz_a2_counts)
            rcp.add_loci_info(_id, nz_het_loci)

        return rcp

    # Generate allele-specific expression effects from read counts.
    def _pr_gen_ase(self, exon_pool=None):
        rcp = self._pr_gen_rcp(exon_pool)

        if len(rcp) != 0:
            aefact = ASEeffectFactory(rcp)
            return aefact.bntest(), aefact.bbtest(), aefact.read_counts

        return None


    def gen_gnm_itvl(self, **kwargs):
        """Generate genomic intervals from Ensembl gene IDs."""
        for _id in self._gene_ids:
            self._pr_get_genome_itvl(_id, **kwargs)

        return self

    def gen_seq(self, shift_factor=0):
        """Generage regulatory sequence matrix."""
        if shift_factor <= 0:
            shift_factor = self.args.shift_factor

        shift = math.ceil(math.sqrt(shift_factor)) ** 2

        if self.mrna_pool:  # Not empty
            self.seq_pool = {
                gene_id: self._pr_gen_ase_seq(itvl=mrna, shift=shift)
                for gene_id, mrna in self.mrna_pool.items()}
        else:
            logging.warning("The mrna_pool is empty, use gen_gnm_itvl() first.")

        return self

    def gen_ase(self):
        """Generate ASE effects."""
        if self.exon_pool:  # Not empty
            self.ase_pool = {
                gene_id: self._pr_gen_ase(exon_pool)
                for gene_id, exon_pool in self.exon_pool.items()}
        else:
            logging.warn("The exon_pool is empty, use gen_gnm_itvl() first.")

        return self

    def save_ase_report(self, opt_file=None):
        """Save the estimated ASE effects into disk."""
        if opt_file is None:
            if self.args.as_ase_report:
                opt_file = self.args.as_ase_report
            else:
                opt_file = self._sample_id + ".ase_report.txt"

        header = "\t".join(["sample_id", "gene_id", "bn_llr", "bn_pval", "bb_llr", "bb_pval",
                            "ase_direction", "snp_info", "snp_phase", "allele_counts"])
        with open(opt_file, "wt") as rfh:
            rfh.write(header + "\n")

            for geneid, ase_effect in self.ase_pool.items():
                if ase_effect:
                    (bnlr, bnpv, _), (bblr, bbpv, dire), (a1_rc, a2_rc, snp_info) = ase_effect

                    _rc, _pos, _phase = ["A1|A2:"], ["CH,PS,ID,REF>ALT:"], ["A1|A2:"]
                    _zip_rec = zip(a1_rc, a2_rc, snp_info)
                    for _a1_rc, _a2_rc, ((chrom, pos, snpid, ref, alt), (a1gt, a2gt)) in _zip_rec:
                        snp_pos = ",".join([str(chrom), str(pos), snpid, ref]) + ">" + alt
                        snp_phase = str(a1gt) + "|" + str(a2gt)
                        snp_rc = str(_a1_rc) + "|" + str(_a2_rc)

                        _rc.append(snp_rc)
                        _pos.append(snp_pos)
                        _phase.append(snp_phase)

                    info_list_1 = [self._sample_id, geneid, bnlr, bnpv, bblr, bbpv, dire]
                    info_list_2 = [x[0] + ";".join(x[1:]) for x in [_pos, _phase, _rc]]
                    est_result = "\t".join([str(x) for x in info_list_1 + info_list_2])
                    rfh.write(est_result + "\n")

        return self

    def save_train_set(self, opt_file=None, save_fmt="fa.gz"):
        """Save seqeunce and ASE effects into disk.

        Args:
            opt_file (str, optional, None): The output file.
            save_fmt (str, optional, fa): output format. ["fa", "fa.gz"]

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
            _ntseq, _ase = self.seq_pool[_gene_id], self.ase_pool[_gene_id]
            _ase_effect = "{}|{}|{}".format(_ase[0][1], *_ase[1][1:]) if _ase else "1|1|0"
            _opt_str += ">{}|{}|{}\n{}\n".format(_sample_id, _gene_id, _ase_effect,
                                                 "\n".join(textwrap.wrap(_ntseq)))

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
    logging.warning("This module should not be executed directly.")
