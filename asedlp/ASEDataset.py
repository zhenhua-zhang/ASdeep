#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# File name : ASEDataset.py
# Author    : Zhenhua Zhang
# E-mail    : zhenhua.zhang217@gmail.com
# Created   : Mon 22 Jun 2020 11:16:49 AM CEST
# Version   : v 0.1.0
# License   : MIT
#
"""A module to load and pre-process sequence matrix."""

import math
import logging

import numpy as np

from torch.utils.data import Dataset, DataLoader, Subset
from statsmodels.sandbox.stats.multicomp import multipletests

from .zutils import logger


class MultipleTestAdjustMent(object):
    """Adjust p-values for multiple tests.

    This class is a dataset-wide transformer.
    """
    def __init__(self, alpha=0.05, method="fdr_bh"):
        self.alpha = 0.05
        self.method = method

    def __call__(self, dataset: list):
        p_val_list = []
        label_list = []
        matrix_list = []

        for matrix, (p_val, label) in dataset:
            p_val_list.append(p_val)
            label_list.append(label)
            matrix_list.append(matrix)

        p_val_list = multipletests(p_val_list, alpha=self.alpha, method=self.method)[1]

        return tuple(zip(matrix_list, tuple(zip(p_val_list, label_list))))


class ReshapeMatrixAndPickupLabel(object):
    """Reshape the sequence matrix into a given size.

    This class is an element-wise transformer.
    Args:
        -1, 0, 1 reprsents ASE effect prone to allele A, no ASE effects, ASE
        Effects prone to allele B in the raw ASE quantification results,
        however, both PyTorch and TensorFlow require labels no less than 0.
    """
    def __init__(self, pthd=0.05):
        self.pthd = pthd

    def __call__(self, sample):
        matrix, labels = sample

        n_channel, length, n_type_nucl = matrix.shape
        length = int(math.sqrt(length * n_type_nucl))

        matrix = matrix.reshape((n_channel, length, length))
        if labels[0] < self.pthd:
            label = labels[1] + 1
        else:
            label = 1

        return (matrix, label)


class ASEDataset(Dataset):
    def __init__(self, gene_id, file_path_pool, element_trans=None, dataset_trans=None):
        """
        Args:
            gene_id   (string): Gene ID (Ensembl gene ID) to train on.
            file_path_pool  (string): Pattern to find the numpy file.
            element_trans (callable, optional): Optional transfrom to be applied
            on a sample.
            dataset_trans (callable, optional): Optional transfrom to be applied
            on the whole dataset.

        NOTE:
            1. In previous implementation, it's not handy to do multiple-test
            adjustment as the sampel was loaded at getitem. Thus, in current
            implementation, the data will be loaded into memory and then further
            operation will be done.
        """
        self.gene_id = gene_id
        self.element_trans = element_trans
        self.dataset_trans = dataset_trans
        self.file_path_pool = file_path_pool

        self.dataset_pool = self._load_data()

    def __len__(self):
        return len(self.file_path_pool)

    def __getitem__(self, idx):
        sample = self.dataset_pool[idx]

        if self.element_trans:
            sample = self.element_trans(sample)

        return sample

    def _load_data(self):
        """Load dataset."""
        temp_list = []
        for idx in range(len(self)):
            gene_id = self.gene_id
            if self.file_path_pool:
                file_path = self.file_path_pool[idx]
            else:
                raise ValueError()
            dataset = np.load(file_path, allow_pickle=True).get(self.gene_id)

            if dataset is None:
                logger.error("[E]: Failed to find {} in file: {}".format(gene_id, file_path))
                return None

            matrix = dataset[0].astype(np.float32)
            length = matrix.shape[1]
            max_sampl_num = int(math.sqrt(length / 2)) ** 2 * 2
            if max_sampl_num != length:
                logger.warn("Reshape the input data for the product of length * width has no integer square root solution")
                matrix = matrix[:, :max_sampl_num, :]

            label = dataset[1]
            sample = (matrix, label)
            temp_list.append(sample)

        if self.dataset_trans:
            temp_list = self.dataset_trans(temp_list)

        return tuple(temp_list)

    def _items(self, idx, labels=True):
        pos = 1 if labels else 0
        if idx is None:
            for idx in range(len(self)):
                yield self[idx][pos]
        else:
            yield self[idx][pos]

    def get_labels(self, idx=None):
        return self._items(idx)

    def get_matrix(self, idx=None):
        return self._items(idx, False)

