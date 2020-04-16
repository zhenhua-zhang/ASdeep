#!/usr/bin/env python3
# -*- utf-8 -*-
import numpy as np

from scipy.stats import beta
from scipy.stats import betabinom
from scipy.stats import chi2
from scipy.optimize import minimize
from scipy.optimize import Bounds

def bblikehood(params, k_vec, n_vec):
    if len(params) != 2:
        raise TypeError("The params should be a list with length 2.")

    bb_a, bb_b = params
    lh_vec = [betabinom.logpmf(k, n, bb_a, bb_b) for k, n in zip(k_vec, n_vec)]

    return -sum(lh_vec)

def llhtest():
    x0 = np.array([0.5, 0.5])
    k_vec = np.array([15, 17, 19, 25, 22, 23])
    n_vec = np.array([39, 43, 48, 58, 50, 45])

    results = minimize(
        bblikehood, x0, (k_vec, n_vec),
        method="BFGS"# , options={"maxfev": 10000}
    )

    est_a, est_b = results["x"]
    estll = results["fun"]

    hyp_a = hyp_b = results["x"].mean()
    hypll = bblikehood((hyp_a, hyp_b), k_vec, n_vec)

    llr = - 2 * (estll - hypll)
    p = chi2.sf(llr, 1)

    print("           Estimated a, b: {}, {}".format(est_a, est_b))
    print("          Hypothysis a, b: {}, {}".format(hyp_a, hyp_b))
    print("Likelihood ratio, p-value: {}, {}".format(llr, p))


import tables

def read_hdf5():
    h5file = tables.open_file(
        "/home/umcg-zzhang/Documents/git/WASP/examples/example_data/haps.h5",
        mode="r", title="haplotype HDF5 table"
    )

    print(h5file)
    h5file.close()


def test():
    import os
    import tables
    import gffutils

    db_path = "chr1_gtf_sql3.db"
    if os.path.exists("./chr1_gtf_sql3.db"):
        db = gffutils.FeatureDB(db_path)
    else:
        fn = "~/Documents/projects/ASECausalSNPPrioritization/inputs/Ensembl_references/Homo_sapiens.GRCh37.75_chr1.gtf"
        db = gffutils.create_db(
            fn, "./chr1_gtf_sql3.db", disable_infer_genes=True,
            disable_infer_transcripts=True
        )

    gene_id = "ENSG00000158828"
    mrna_pool = db.children(gene_id, featuretype="transcript")

    mrna = next(mrna_pool, False)
    if mrna:
        parent_id = mrna.id
    else:
        parent_id = gene_id

    exon_pool = db.children(parent_id, featuretype="exon")

    ref_read_file = "/home/umcg-zzhang/Documents/projects/ASECausalSNPPrioritization/workdir/optdir/AC47H5ACXX-3-16/waspOptdir/perChrom/1/AC47H5ACXX-3-16_1.refAlleleCounts.h5"
    ref_read_pool = tables.open_file(ref_read_file)
    ref_read_chr1 = ref_read_pool.get_node("/1")

    alt_read_file = "/home/umcg-zzhang/Documents/projects/ASECausalSNPPrioritization/workdir/optdir/AC47H5ACXX-3-16/waspOptdir/perChrom/1/AC47H5ACXX-3-16_1.altAlleleCounts.h5"
    alt_read_pool = tables.open_file(alt_read_file)
    alt_read_chr1 = alt_read_pool.get_node("/1")

    other_read_file = "/home/umcg-zzhang/Documents/projects/ASECausalSNPPrioritization/workdir/optdir/AC47H5ACXX-3-16/waspOptdir/perChrom/1/AC47H5ACXX-3-16_1.otherAlleleCounts.h5"
    other_read_pool = tables.open_file(other_read_file)
    other_read_chr1 = other_read_pool.get_node("/1")

    if mrna:
        print("mRNA coord: {}-{}".format(mrna.start, mrna.stop))
        ref_read_counts = ref_read_chr1[mrna.start: mrna.stop] 
        alt_read_counts = alt_read_chr1[mrna.start: mrna.stop] 
    else:
        print("it looks there's no genomic interval for the gene ID.")

    for id, exon in enumerate(exon_pool):
        print(exon.astuple())
        print("exon coord: {}, {}:{}-{}".format(exon.id, exon.seqid, exon.start, exon.stop))
        ref_read_counts = ref_read_chr1[exon.start: exon.stop]
        alt_read_counts = alt_read_chr1[exon.start: exon.stop]
        other_read_counts = other_read_chr1[exon.start: exon.stop]

        print(ref_read_counts[ref_read_counts > 0])
        print(alt_read_counts[alt_read_counts > 0])
        print(other_read_counts[other_read_counts > 0])

    seq_tab = tables.open_file("./../../workdir/fastah5db/sequence.h5")
    hap_tab = tables.open_file("./../../workdir/snph5db/1/haplotype.h5")
    snp_tab = tables.open_file("./../../workdir/snph5db/1/snps_tab.h5")
    snp_idx = tables.open_file("./../../workdir/snph5db/1/snps_index.h5")

    seq_tab_chr1 = seq_tab.get_node("/1")
    snp_idx_chr1 = snp_idx.get_node("/1")

    upstream = seq_tab_chr1[mrna.start - 1000: mrna.start]
    upstream_snps = snp_idx_chr1[mrna.start - 1000: mrna.start]

    hap_chr1 = hap_tab.get_node("/1")
    hap_chr1_phase_1 = hap_tab.get_node("/phase_1")

    hap_samples = hap_tab.get_node("/samples_1")
    sampleid_to_idx = {sample[0].decode(): idx for idx, sample in enumerate(hap_samples)}

    snp_tab_chr1 = snp_tab.get_node("/1")

    for _a1, _a2 in list(zip(upstream, upstream_snps)):
        if _a2 == -1:
            pass
        else:
            snp_info = snp_tab_chr1[_a2]
            hap_info = hap_chr1[_a2, 0:2]
            a1_idx, a2_idx = hap_info
            a1, a2 = snp_info[a1_idx+2], snp_info[a2_idx+2]

    seq_tab.close()
    hap_tab.close()
    snp_tab.close()
    snp_idx.close()
    ref_read_pool.close()
    alt_read_pool.close()
    other_read_pool.close()


from collections import UserDict
import itertools

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

    def unpack(self):
        ref_read_pool, alt_read_pool = [], []
        for _key, _value in self.data.items():
            ref_read_pool.extend(_value[0])
            alt_read_pool.extend(_value[1])
        return ref_read_pool, alt_read_pool


def flatten(iptl):
    out_list = []
    for ele in iptl:
        if isinstance(ele, (list, tuple)):
            out_list.extend(flatten(ele))
        else:
            out_list.append(ele)
    return out_list

