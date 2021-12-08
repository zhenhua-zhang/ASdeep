# Author : Zhenhua Zhang
# E-mail : zhenhua.zhang217@gmail.com
# Created: Mar 12, 2020
# Updated: Oct 06, 2021

import numpy as np
import h5py as h5
from torch.utils.data import Dataset

from .database import HilbertCurve
from .zutils import LogManager


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
        if not isinstance(record, tuple) or len(record) != 2:
            raise TypeError("The input should be a tuple of length 2.")

        hbcurve, attrs = record

        if self._n_bp > 0:
            strand = attrs.get(self._strand_key, 1)
            hbcurve = HilbertCurve(hbcurve, self._kmers, strand)
            hbcurve = hbcurve.subset(self._n_bp)

        return hbcurve, attrs

class MaskHomoSites:
    """Mask homozygous sites in the sequnece."""
    def __init__(self, flank: int = 25, fills = -1, label_key: str = "ASE"):
        self._flank = flank
        self._fills = fills
        self._label_key = label_key

    def __call__(self, record: tuple):
        if not isinstance(record, tuple) or len(record) != 2:
            raise TypeError("The input should be a tuple of length 2.")

        hbcurve, attrs = record

        if isinstance(hbcurve, np.ndarray) or isinstance(hbcurve, str):
            hbcurve = HilbertCurve(hbcurve)

        hbcurve.mask_homo(self._flank, self._fills)
        return hbcurve.hbcmat, attrs


class XyTransformer:
    """Split each h5.Dataset into X (matrix) and y (label)."""
    def __init__(self, adj_label=lambda x: int(x) + 1, label_key: str = "ASE"):
        self._adj_label = adj_label
        self._label_key = label_key

    def __call__(self, record: tuple):
        if not isinstance(record, tuple) or len(record) != 2:
            raise TypeError("The input should be a tuple of length 2.")

        hbcurve, attrs = record
        label = attrs.get(self._label_key, None)
        if label is not None:
            label = self._adj_label(label)

        if isinstance(hbcurve, HilbertCurve):
            hbcurve = hbcurve.hbcmat

        return hbcurve, label


class ASEDataset(Dataset):
    """ASE dataset.

    Attributes:
        samples: All availabel samples in the database.

    Notes:
        1. ASEDataset inherites PyTorch's Dataset which accelerates data
        loading by multiprocessing, however, HDF5 database does not support
        multiprocessing data loading schemes after opening the file. The
        solution is to open the file for every loading, which could introduce
        extra resource burden.
    """
    def __init__(self, dbpath: str, transformers=None, n_cpus: int = 1,
                 label_key: str = "ASE",
                 logman: LogManager = LogManager("ASEDataset")):
        self._n_cpus = n_cpus
        self._dbpath = dbpath
        self._transformers = transformers
        self._label_key = label_key
        self._logman = logman

        with h5.File(dbpath, "r") as database:
            self._samples = list(database.keys())

    def __len__(self):
        return len(self._samples)

    def __getitem__(self, idx):
        if isinstance(idx, int):
            if idx < 0:
                raise IndexError("Index should be >= 0.")

            idx = f"/{self._samples[idx]}"
        elif isinstance(idx, str):
            idx = f"/{idx.strip('/')}"

        self._logman.debug(idx)
        with h5.File(self._dbpath, "r") as database:
            record = database.get(idx)
            record = (record[()], dict(record.attrs))
            _hbcurve, _attrs = self._transform(record, self._transformers)

        return (_hbcurve, _attrs)

    def __contains__(self, key):
        return key in self._samples

    @staticmethod
    def _transform(record, transformers):
        if not isinstance(transformers, (list, tuple)):
            transformers = [transformers]

        for per_trans in transformers:
            record = per_trans(record)

        return record

    def _items(self, idx=None, labels=True, apply_trans=None):
        pos = 1 if labels else 0

        if idx is None:
            idx = self._samples
        elif isinstance(idx, int) and 0 <= idx < len(self):
            idx = [self._samples[idx]]
        elif isinstance(idx, str) and idx in self:
            idx = [idx]
        else:
            raise KeyError()

        with h5.File(self._dbpath, "r") as database:
            for per_idx in idx:
                record = database.get(per_idx)
                record = (record[()], dict(record.attrs))
                if apply_trans:
                    yield self._transform(record, apply_trans)[pos]
                else:
                    yield self._transform(record, self._transformers)[pos]

    @property
    def samples(self):
        return self._samples

    def get_labels(self, idx=None):
        """Get labels."""
        return self._items(idx, apply_trans=XyTransformer())

    def get_matrix(self, idx=None):
        """Get matrix."""
        return self._items(idx, labels=False)

