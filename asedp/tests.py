#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import unittest

import numpy as np

try:
    from ASEFactory import ASEFactory
    from ASEDataset import ASEDataset, SeqToHelbertAndMakeLabel
    from HelbertCurve import HelbertCurve
except ImportError as ime:
    print("ImportError: {}".format(ime), file=sys.stderr)
    sys.exit()


class Args:
    pass

class TestASEFactory(unittest.TestCase):

    args = Args()
    # Simulate commandline options.
    setattr(args, 'variants_file', '../data/test.vcf.gz')
    setattr(args, 'genome_seq_file', '../data/test.fa')
    setattr(args, 'gene_feature_file', '../data/test.gtf')
    setattr(args, 'ase_readcount_file', '../data/test_wasp_ase_readcounts.ssv')
    setattr(args, 'gene_ids', ['ENSG00000198062'])
    setattr(args, 'sample_id', 'individual_1')
    setattr(args, 'gene_id_file', None)
    setattr(args, 'as_ase_report', '../data/test_ase_report.txt')
    setattr(args, 'as_train_set', '../data/test_ase_train_set.fa.gz')

    def test_ASEFactory(self):
        ASEFactory(self.args) \
                .gen_gnm_itvl() \
                .gen_seq() \
                .gen_ase() \
                .save_ase_report() \
                .save_train_set()


class TestASEDataset(unittest.TestCase):
    gene_id = "ENSG00000188976"
    file_pat = ["../examples/chr1.ntsq_and_ase.fa.gz"]

    def test_ASEDataset(self):
        element_trans = SeqToHelbertAndMakeLabel(show_pic=True)
        dataset = ASEDataset(self.gene_id, self.file_pat,
                             element_trans=element_trans)
        self.assertEqual(len(dataset), len(dataset.file_path_pool))
        self.assertEqual(dataset.gene_id, self.gene_id)
        self.assertEqual(dataset[0][1], 1)


class TestHelbertCurve(unittest.TestCase):
    kmers = 4

    def test_helbert_curve_all_base(self):
        seq = "ATCGMRSYKWmrsykwATCGMRSYKWmrsykw"
        ohhc = HelbertCurve(seq, kmers=self.kmers) \
                .seq_to_hcurve() \
                .hcurve_to_img(output_prefix="all_base.", figsize=10) \
                .get_onehot_hcurve()

        self.assertEqual(np.sum(ohhc), 58)

    def test_helbert_curve_htz_base(self):
        seq = "ATCGWCTATGACGTAAATCGGCTATGACGTAA"
        ohhc = HelbertCurve(seq, kmers=self.kmers) \
                .seq_to_hcurve() \
                .hcurve_to_img(output_prefix="het_base.", figsize=10) \
                .get_onehot_hcurve()

        self.assertEqual(np.sum(ohhc), 58)

    def test_helbert_curve_sym_base(self):
        seq = "ATCGGCTATGACGTAAATCGGCTATGACGTAA"
        ohhc = HelbertCurve(seq, kmers=self.kmers) \
                .seq_to_hcurve() \
                .hcurve_to_img(output_prefix="sym_base.", figsize=10) \
                .get_onehot_hcurve()

        self.assertEqual(np.sum(ohhc), 58)

    def test_helbert_curve(self):
        seq = "ACGATGCTCAG"
        ohhc = HelbertCurve(seq, kmers=self.kmers) \
                .seq_to_hcurve() \
                .hcurve_to_img(output_prefix="sixteen_base.", figsize=2.5) \
                .get_onehot_hcurve()
        self.assertGreater(np.sum(ohhc), 0)


if __name__ == "__main__":
    unittest.main()
