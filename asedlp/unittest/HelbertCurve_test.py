#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author      : Zhenhua Zhang
# Email       : zhenhua.zhang217@gmail.com
# License     : MIT
# Created date: 2020 Jul 17 14:25:15

"""Unit test for `HelbertCurve` class in asedlp.HelbertCurve module"""

import sys
import unittest
from textwrap import wrap

import numpy as np
sys.path.insert(1, "..")

try:
    from HelbertCurve import HelbertCurve
    from zutils import logger
except ImportError as err:
    logger.error(err)
    sys.exit()


class TestHelbertCurve(unittest.TestCase):
    kmers = 4

    def test_helbert_curve_all_base(self):
        seq = "ATCGMRSYKWmrsykwATCGMRSYKWmrsykw"
        ohhc = HelbertCurve(seq, kmers=self.kmers) \
                .seq_to_hcurve() \
                .hcurve_to_img(output_prefix="all_base.", figsize=10) \
                .get_onehot_hcurve()

        self.assertEqual(np.sum(ohhc), 58)

        _seq = "\n".join(wrap(seq, 80))
        # logger.info("Found {} {}-mers:\n{}".format(np.sum(ohhc), self.kmers, _seq))

    def test_helbert_curve_htz_base(self):
        seq = "ATCGWCTATGACGTAAATCGGCTATGACGTAA"
        ohhc = HelbertCurve(seq, kmers=self.kmers) \
                .seq_to_hcurve() \
                .hcurve_to_img(output_prefix="het_base.", figsize=10) \
                .get_onehot_hcurve()

        self.assertEqual(np.sum(ohhc), 58)

        _seq = "\n".join(wrap(seq, 80))
        # logger.info("Found {} {}-mers:\n{}".format(np.sum(ohhc), self.kmers, _seq))

    def test_helbert_curve_sym_base(self):
        seq = "ATCGGCTATGACGTAAATCGGCTATGACGTAA"
        ohhc = HelbertCurve(seq, kmers=self.kmers) \
                .seq_to_hcurve() \
                .hcurve_to_img(output_prefix="sym_base.", figsize=10) \
                .get_onehot_hcurve()

        self.assertEqual(np.sum(ohhc), 58)

        _seq = "\n".join(wrap(seq, 80))
        # logger.info("Found {} {}-mers:\n{}".format(np.sum(ohhc), self.kmers, _seq))


if __name__ == "__main__":
    unittest.main()
