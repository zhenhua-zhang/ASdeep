#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# File name: ASEDataset.py
# Author   : Zhenhua Zhang
# E-mail   : zhenhua.zhang217@gmail.com
# Created  : Mon 22 Jun 2020 11:16:49 AM CEST
# License  : MIT
#
"""A module to load and pre-process sequence matrix."""

import math
import logging
from multiprocessing import Pool

import pyfaidx
from torch.utils.data import Dataset

from .hbcurve import HilbertCurve
from .zutils import fdr_bh


class MultipleTestAdjustment:
    """Adjust p-values for multiple tests.

    This class is a dataset-wide transformer.
    """
    def __init__(self, alpha: float = 0.05, method: str = "fdr_bh"):
        self.alpha = alpha
        self.method = method

    def __call__(self, dataset: list):
        pval_pool = []
        mtrx_pool = []
        label_pool = []

        for item in dataset:
            if len(item) == 1:
                (matrix, (p_val, label)) = item[0]
                pval_pool.append(p_val)
                mtrx_pool.append(matrix)
                label_pool.append(label)
            else:
                logging.warning("The length of input dataset should be in length 1.")

        pval_pool = fdr_bh(pval_pool)

        return tuple(zip(mtrx_pool, tuple(zip(pval_pool, label_pool))))


class ReshapeMatrixAndPickupLabel:
    """Reshape the sequence matrix into a given size.

    This class is an element-wise transformer.

    Attributes:
        -1, 0, 1 reprsents ASE effect prone to allele A, no ASE effects, ASE
        Effects prone to allele B in the raw ASE quantification results,
        however, both PyTorch and TensorFlow require labels no less than 0.
    """
    def __init__(self, pthd=0.05):
        self.pthd = pthd

    def __call__(self, sample):
        _matrix, labels = sample

        n_channel, length, n_type_nucl = _matrix.shape
        length = int(math.sqrt(length * n_type_nucl))

        _matrix = _matrix.reshape((n_channel, length, length))
        _label = labels[1] + 1 if labels[0] < self.pthd else 1

        return (_matrix, _label)


class SeqToHilbertAndMakeLabel:
    """Convert sequence into matrix of Hilbert curve and fit label to training.

    This class is a dataset-wide transformer.
    """
    def __init__(self, pthd=0.05, show_pic=False, onehot=False, matrix=True):
        self.pthd = pthd
        self.show_pic = show_pic
        self.onehot = onehot
        self.matrix = matrix

    def __call__(self, sample, **kwargs):
        sequence, (p_val, label) = sample

        if p_val is None or label is None:
            _p_val = 1
            _label = -1
        else:
            _p_val = p_val
            if p_val <= self.pthd:
                _label = label + 1
            else:
                _label = 1

        _hcurve = HilbertCurve(sequence, **kwargs).seq_to_hcurve()

        if self.show_pic:
            _hcurve = _hcurve.hcurve_to_img()

        if self.matrix:
            _matrix = _hcurve.get_hcurve(self.onehot)
            return (_matrix.astype(float), (_p_val, _label))

        return (_hcurve, (_p_val, _label))


class ASEDataset(Dataset):
    """ASE dataset.

    Attributes:
        gene_id (string): Gene ID (Ensembl gene ID) to train on.
        file_path_pool (string): Pattern to find the numpy file.
        element_trans (callable, optional): Optional transfrom to be applied on a sample.
        dataset_trans (callable, optional): Optional transform to be applied on the whole dataset.
        use_bb_pval (bool, optional): Use Beta-Binomial p-value. Default: False

    NOTE:
        1. In previous implementation, it's not handy to do multiple-test
        adjustment as the sample was loaded at __getitem__(). Thus, in current
        implementation, the data will be loaded into memory and then further
        operation will be done.
    """
    def __init__(self,
                 gene_id,
                 file_path_pool,
                 element_trans=None,
                 dataset_trans=None,
                 use_bb_pval=False,
                 num_workers=1):
        self.gene_id = gene_id

        if isinstance(gene_id, str):
            self.gene_id = [gene_id]
        elif not isinstance(gene_id, (list, tuple)):
            raise ValueError("gene_id should be a list/tuple/str!")

        self.element_trans = element_trans
        self.dataset_trans = dataset_trans
        self.file_path_pool = file_path_pool
        self.use_bb_pval = use_bb_pval

        self.num_workers = num_workers
        self.dataset_pool = []

        self._load_data()

    def __len__(self):
        return len(self.dataset_pool)

    def __getitem__(self, idx):
        return self.dataset_pool[idx]

    def _load_data_base(self, file_path):
        seq_pool = pyfaidx.Fasta(file_path)
        chosen_records = [
            seq_pool[idx] for idx in seq_pool.keys()
            if any([_gene_id in idx for _gene_id in self.gene_id])
        ]

        tmp_pool = []
        for record in chosen_records:
            if record is None or len(record) == 0:
                _err_msg = "No \"{}\" in \"{}\"".format(
                    self.gene_id, file_path)
                logging.warning(_err_msg)
                record = (None, (None, None))
            else:
                record_name_list = record.name.split("|")
                if len(record_name_list) > 2:
                    p_val_bn, p_val_bb, label = record_name_list[2:5]
                    p_val = p_val_bb if self.use_bb_pval else p_val_bn
                    p_val, label = float(p_val), int(label)
                else:
                    p_val, label = None, None
                record = (str(record[:].seq), (p_val, label))

            if self.element_trans:
                record = self.element_trans(record)
            tmp_pool.append(record)
        seq_pool.close()

        return tmp_pool

    def _load_data(self):
        with Pool(self.num_workers) as proc_pool:
            self.dataset_pool = proc_pool.map(self._load_data_base, self.file_path_pool)

        if self.dataset_trans:
            self.dataset_pool = self.dataset_trans(self.dataset_pool)

    def _items(self, idx=None, labels=True):
        # Yield items.
        pos = 1 if labels else 0
        if idx is None:
            for _idx, _ in enumerate(self):
                yield self[_idx][pos][1]
        else:
            yield self[idx][pos][1]

    def get_labels(self, idx=None):
        """Get labels."""
        return self._items(idx)

    def get_matrix(self, idx=None):
        """Get matrix."""
        return self._items(idx, False)


if __name__ == "__main__":
    logging.warning("This module should not be executed directly.")