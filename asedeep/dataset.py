#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# File name: ASEDataset.py
# Author   : Zhenhua Zhang
# E-mail   : zhenhua.zhang217@gmail.com
# Created  : Mon 22 Jun 2020 11:16:49 AM CEST
# License  : MIT
#
'''A module to load and pre-process sequence matrix.'''

import math
import logging

import pyfaidx
from torch.utils.data import Dataset
from statsmodels.sandbox.stats.multicomp import multipletests

from .hbcurve import HelbertCurve


class MultipleTestAdjustment:
    '''Adjust p-values for multiple tests.

    This class is a dataset-wide transformer.
    '''

    def __init__(self, alpha: float = 0.05, method: str = 'fdr_bh'):
        self.alpha = alpha
        self.method = method

    def __call__(self, dataset: list):
        pval_pool = []
        mtrx_pool = []
        label_pool = []

        for matrix, (p_val, label) in dataset:
            pval_pool.append(p_val)
            mtrx_pool.append(matrix)
            label_pool.append(label)

        pval_pool = multipletests(pval_pool, alpha=self.alpha, method=self.method)

        return tuple(zip(mtrx_pool, tuple(zip(pval_pool[1], label_pool))))


class ReshapeMatrixAndPickupLabel:
    '''Reshape the sequence matrix into a given size.

    This class is an element-wise transformer.
    Args:
        -1, 0, 1 reprsents ASE effect prone to allele A, no ASE effects, ASE
        Effects prone to allele B in the raw ASE quantification results,
        however, both PyTorch and TensorFlow require labels no less than 0.
    '''

    def __init__(self, pthd=0.05):
        self.pthd = pthd

    def __call__(self, sample):
        _matrix, labels = sample

        n_channel, length, n_type_nucl = _matrix.shape
        length = int(math.sqrt(length * n_type_nucl))

        _matrix = _matrix.reshape((n_channel, length, length))
        _label = labels[1] + 1 if labels[0] < self.pthd else 1

        return (_matrix, _label)


class SeqToHelbertAndMakeLabel:
    '''Convert sequence into matrix of Helbert curve and fit label to training.

    This class is a dataset-wide transformer.
    '''
    def __init__(self, pthd=0.05, show_pic=False, onehot=False, matrix=True):
        self.pthd = pthd
        self.show_pic = show_pic
        self.onehot = onehot
        self.matrix = matrix

    def __call__(self, sample, **kwargs):
        sequence, labels = sample

        if labels[0] <= self.pthd:
            _label = labels[1] + 1
        else:
            _label = 1

        _hcurve = HelbertCurve(sequence, **kwargs).seq_to_hcurve()

        if self.show_pic:
            _hcurve = _hcurve.hcurve_to_img()

        if self.matrix:
            _matrix = _hcurve.get_hcurve(self.onehot)
            return (_matrix.astype(float), _label)

        return (_hcurve, _label)


class ASEDataset(Dataset):
    '''ASE dataset.
    Args:
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
    '''
    def __init__(self, gene_id, file_path_pool, element_trans=None, dataset_trans=None,
                 use_bb_pval=False):
        self.gene_id = gene_id
        self.element_trans = element_trans
        self.dataset_trans = dataset_trans
        self.file_path_pool = file_path_pool
        self.use_bb_pval = use_bb_pval

        self.dataset_pool = self._load_data()

    def __len__(self):
        return len(self.file_path_pool)

    def __getitem__(self, idx):
        return self.element_trans(self.dataset_pool[idx]) if self.element_trans \
            else self.dataset_pool[idx]

    def _load_data(self):
        # Load dataset.
        temp_list = []
        for idx in range(len(self)):
            file_path = self.file_path_pool[idx]

            seq_pool = pyfaidx.Fasta(file_path)
            record = [seq_pool[record_id] for record_id in seq_pool.keys()
                      if self.gene_id in record_id][0]

            if record is None or len(record) == 0:
                _err_msg = 'No \'{}\' in \'{}\''.format(self.gene_id, file_path)
                logging.warning(_err_msg)
                temp_list.append([None, None])
            else:
                p_val_bn, p_val_bb, label = record.name.split('|')[2:5]
                if self.use_bb_pval:
                    p_val = p_val_bb
                else:
                    p_val = p_val_bn

                temp_list.append((str(record[0:].seq), (float(p_val), int(label))))

        if self.dataset_trans:
            return tuple(self.dataset_trans(temp_list))

        return tuple(temp_list)

    def _items(self, idx=None, labels=True):
        # Yield items.
        pos = 1 if labels else 0
        if idx is None:
            for _idx, _ in enumerate(self):
                yield self[_idx][pos]
        else:
            yield self[idx][pos]

    def get_labels(self, idx=None):
        '''Get labels.'''
        return self._items(idx)

    def get_matrix(self, idx=None):
        '''Get matrix.'''
        return self._items(idx, False)


if __name__ == "__main__":
    logging.warning("This module should not be executed directly.")
