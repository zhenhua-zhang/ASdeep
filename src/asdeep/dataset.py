# Author : Zhenhua Zhang
# E-mail : zhenhua.zhang217@gmail.com
# Created: Mar 12, 2020
# Updated: Oct 06, 2021

import math

import h5py as h5
import numpy as np
from hilbert import decode
from torch.utils.data import Dataset

from .database import HilbertCurve
from .zutils import LogManager
from .zutils import fdr_bh


class MultipleTestAdjustment:
    """Adjust p-values for multiple tests.

    This class is a dataset-wide transformer.
    """
    def __init__(self,
                 alpha: float = 0.05,
                 method: str = "fdr_bh",
                 logman: LogManager = LogManager("MultipleTestAdjustment")):
        self.alpha = alpha
        self.method = method
        self.logman = logman

    def __call__(self, dataset: list):
        pval_pool, mtrx_pool, label_pool = [], [], []

        for item in dataset:
            if len(item) == 1:
                (matrix, (p_val, label)) = item[0]
                pval_pool.append(p_val)
                mtrx_pool.append(matrix)
                label_pool.append(label)
            else:
                self.logman.warning("The length of input dataset should be 1.")

        if self.method == "fdr_bh":
            pval_pool = fdr_bh(pval_pool)
        else:
            raise ValueError("Uknown mutliple test adjustment method: " +
                             "{}".format(self.method))

        return tuple(zip(mtrx_pool, tuple(zip(pval_pool, label_pool))))


class ReshapeMatrixAndPickupLabel:
    """Reshape the sequence matrix into a given size.

    This class is an element-wise transformer.

    Notes:
        -1, 0, 1 reprsents ASE effect prone to allele A, no ASE effects, ASE
        Effects prone to allele B in the raw ASE quantification results,
        however, both PyTorch and TensorFlow require labels no less than 0.
    """
    def __init__(self, pthd: float = 0.05):
        self.pthd = pthd

    def __call__(self, sample):
        hbcurve, labels = sample

        matrix = hbcurve.hbcurve
        n_channel, length, n_type_nucl = matrix.shape
        length = int(math.sqrt(length * n_type_nucl))

        matrix = matrix.reshape((n_channel, length, length))
        label = labels[1] + 1 if labels[0] < self.pthd else 1

        return (matrix, label)


class SubsetHilbertCurve:
    """Subset a Hilbert curve

    NOTE: the kmers should be made using the same parameters
    """
    def __init__(self, n_bp: int = 0, strand_key: str = "strand",
                 kmers: int = 4, pad_by=-1,
                 logman: LogManager = LogManager("SubsetHilbertCurve")):
        self._n_bp = n_bp
        self._strand_key = strand_key
        self._kmers = kmers
        self._pad_by = pad_by
        self._logman = logman

    def __call__(self, record: tuple):
        hbcurve, attrs = record
        strand = attrs.get(self._strand_key, 1)

        if self._n_bp > 0:
            hbcurve = HilbertCurve(hbcurve).subset(self._n_bp, strand)

        return hbcurve, attrs

class XyTransformer:
    """Split each h5.Dataset into X (matrix) and y (label)."""
    def __init__(self, adj_label=lambda x: int(x) + 1, label_key: str = "ASE"):
        self._adj_label = adj_label
        self._label_key = label_key

    def __call__(self, record: tuple):
        hbcurve, attrs = record
        label = attrs.get(self._label_key, None)
        if label is not None:
            label = self._adj_label(label)
        return hbcurve.hbcurve, label


class ASEDataset(Dataset):
    """ASE dataset.

    Attributes:
    """
    def __init__(self, database: h5.File, n_cpus: int = 1, transformers=None,
                 label_key: str = "ASE",
                 logman: LogManager = LogManager("ASEDataset")):
        self._n_cpus = n_cpus
        self._database = database
        self._transformers = transformers
        self._label_key = label_key
        self._logman = logman

        self._samples = list(self._database.keys())

    def __len__(self):
        return len(self._samples)

    def __getitem__(self, idx):
        if isinstance(idx, int):
            if idx < 0: raise IndexError("Index should be >= 0.")
            idx = f"/{self._samples[idx]}"
        elif isinstance(idx, str):
            if not idx.startswith("/"):
                idx = f"/{idx}"

        record = self._database.get(idx)
        record = (record[()], record.attrs)
        _hbcurve, _attrs = self._apply_transformer(record, self._transformers)
        return (_hbcurve, _attrs)

    def __contains__(self, key):
        return key in self._samples

    @staticmethod
    def _apply_transformer(record, transformers):
        if not isinstance(transformers, (list, tuple)):
            transformers = [transformers]

        for per_trans in transformers:
            record = per_trans(record)

        return record

    def _items(self, idx=None, labels=True): # Yield items.
        pos = 1 if labels else 0
        if idx is None:
            for idx_ in range(len(self)):
                yield self[idx_][pos]
        else:
            yield self[idx][pos]

    @property
    def samples(self):
        return self._samples

    def get_labels(self, idx=None):
        """Get labels."""
        return self._items(idx)

    def get_matrix(self, idx=None):
        """Get matrix."""
        return self._items(idx, False)

