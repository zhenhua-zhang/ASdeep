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
    from asedlp.ASEDataset import ASEDataset
    from asedlp.DLPFactory import DLPFactory
    from asedlp.DLPFactory import CNNModel
except ImportError as ime:
    print("ImportError: {}".format(ime), file=sys.stderr)


class TestASEDLP(unittest.TestCase):
    gene_id = "ENSG00000185220"
    file_pat = "../examples/**/*_1_matrix_and_ase.npz"
    def test_ASEDataset(self):
        dataset = ASEDataset(self.gene_id, self.file_pat)
        self.assertEqual(len(dataset), len(dataset.file_path_pool))
        self.assertEqual(dataset.gene_id, self.gene_id)

    def test_DLPFactory(self):
        net = CNNModel()
        factory = DLPFactory(net, self.gene_id, self.file_pat)
        factory.init() \
                .load_dataset(batch_size=3) \
                .train()


if __name__ == "__main__":
    unittest.main()
