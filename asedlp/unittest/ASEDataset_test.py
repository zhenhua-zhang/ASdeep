#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author      : Zhenhua Zhang
# Email       : zhenhua.zhang217@gmail.com
# License     : MIT
# Created date: Sun 22 Mar 2020 11:18:52 AM CET
# Last update : Sun 22 Mar 2020 11:19:12 AM CET

"""Unit test for `ASEDataset` class in `asedlp.ASEDataset` module"""

import sys
import unittest

sys.path.append("..")

try:
    from ASEDataset import ASEDataset, SeqToHelbertAndMakeLabel
except ImportError as ime:
    print("ImportError: {}".format(ime), file=sys.stderr)


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

# class TestDLPFactory(unittest.TestCase):
    # def test_DLPFactory(self):
    #     net = CNNModel()
    #     factory = DLPFactory(net, self.gene_id, self.file_pat)
    #     factory.init() \
    #             .load_dataset(batch_size=3) \
    #             .train()


if __name__ == "__main__":
    unittest.main()
