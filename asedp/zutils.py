#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# File name : zutils.py
# Author        : Zhenhua Zhang
# E-mail        : zhenhua.zhang217@gmail.com
# Created     : Mon 22 Jun 2020 11:25:03 AM CEST
# Version     : v0.1.0
# License     : MIT
#

import sys
import copy
import numpy as np

from collections import UserDict

# Reference: https://genomevolution.org/wiki/index.php/Ambiguous_nucleotide
M2B = { # monoallelic base-pair to biallelic base.
    "AA": "A", "AC": "M", "AG": "R", "AT": "W",
    "CA": "m", "CC": "C", "CG": "S", "CT": "Y",
    "GA": "r", "GC": "s", "GG": "G", "GT": "K",
    "TA": "w", "TC": "y", "TG": "k", "TT": "T",
    "NN": "N"
}

# M2B = { # monoallelic base-pair to biallelic base.
#     (65, 65): "A", (65, 67): "M", (65, 71): "R", (65, 84): "W",
#     (67, 65): "m", (67, 67): "C", (67, 71): "S", (67, 84): "Y",
#     (71, 65): "r", (71, 67): "s", (71, 71): "G", (71, 84): "K",
#     (84, 65): "w", (84, 67): "y", (84, 71): "k", (84, 84): "T",
#     (78, 78): "N",
# }

B2M = { # Biallelic base to monoallelic base-pair
    "A": "AA", "M": "AC", "R": "AG", "W": "AT",
    "m": "CA", "C": "CC", "S": "CG", "Y": "CT",
    "r": "GA", "s": "GC", "G": "GG", "K": "GT",
    "w": "TA", "y": "TC", "k": "TG", "T": "TT",
    "N": "NN"
}


def print_arguments(arguments, fwidth=None, to_head=("subcommand",)):
    """Print artuments from command lines"""
    print("Arguments for current run: ", file=sys.stderr)
    _dict = vars(arguments)
    _pair = [
        (_dst, _arg) for _dst, _arg in _dict.items() if _dst not in to_head
    ]
    _pair = sorted(_pair, key=lambda x: len(x[0]))

    _to_head = [
        (_dst, _arg) for _dst, _arg in _dict.items() if _dst in to_head
    ]

    for prior in _to_head:
        _pair.insert(0, prior)

    if fwidth is None:
        fwidth = len(_pair[-1][0]) + 1

    for _dst, _arg in _pair:
        if _arg is None:
            _arg = "None"

        if isinstance(_arg, (list, tuple)):
            _arg = ", ".join(map(str, _arg))

        print(
            "    {d: <{w}}: {a: <{w}}".format(d=_dst, a=_arg, w=fwidth),
            file=sys.stderr
        )


class SmartDict(UserDict):
    """Container for the data set.
    """

    def __init__(self):
        super().__init__()
        self.data = {}

    def __setitem__(self, key, val):
        self.data[key] = val

    def __getitem__(self, key):
        return self.data[key]

    def __repr__(self):
        return repr(self.__dict__)

    def __delitem__(self, key):
        del self.data[key]

    def clear(self):
        return self.data.clear()

    def copy(self):
        return self.data.copy()


def cmp(a, b):
    return (a > b) - (a < b)


def flatten(l):
    """Flatten a list recursively.
    Args:
        l (list): The list to be flatten.
    Returns:
        out_list (list): The flatten list.
    """
    out_list = []
    for x in l:
        if isinstance(x, (list, tuple)):
            out_list.extend(flatten(x))
        else:
            out_list.append(x)
    return out_list


def insert_or_append(dict1, dictlike2):
    dict1 = copy.deepcopy(dict1)

    if not hasattr(dictlike2, "__getitem__") or not hasattr(dictlike2, "items"):
        raise AttributeError(
            "dictlike2 should at least has "
            "`__getitem__()`, aka `[]`, methods and itesm()")

    for key, val in dictlike2.items():
        if key in dict1:
            dict1[key].append(val)
        else:
            dict1[key] = [val]

    return dict1

# def split_train_test(dataset: tf.data.Dataset, test_prop: float = 0.3,
#         shuffle: bool = True, shuf_buf_size: int = 32,
#         batch_size: int = 16) -> (tf.data.Dataset, tf.data.Dataset):
#     """Split a dataset into train and test dataset.
# 
#     Source: https://stackoverflow.com/questions/51125266/how-do-i-split-tensorflow-datasets
#     """
#     border = int(10 * test_prop)
# 
#     recover = lambda x, y: y
#     is_test = lambda x, y: x % 10 < border
#     is_train = lambda x, y: x % 10 >= border
# 
#     if shuffle:
#         dataset = dataset.shuffle(shuf_buf_size)
#     train_dtst = dataset.enumerate() \
#             .filter(is_train) \
#             .map(recover) \
#             .batch(batch_size) \
#             .cache()
# 
#     test_dtst = dataset.enumerate() \
#             .filter(is_test) \
#             .map(recover) \
#             .batch(batch_size) \
#             .cache()
# 
#     return train_dtst, test_dtst

def make_gif(fp_in, fp_out, duration=400, loop=0):
    import glob
    from PIL import Image

    try:
        fp_in_imgs = glob.glob(fp_in)
        img, *imgs = [Image.open(f) for f in sorted(fp_in_imgs)]
        img.save(fp=fp_out, format="GIF", append_images=imgs, save_all=True,
                duration=duration, loop=loop)
    except ValueError as err:
        print(err)


# this code is copied from statsmodel source code, to avoid installation of statsmodel for computing fdr-corrected p-values
# link: http://statsmodels.sourceforge.net/ipdirective/_modules/scikits/statsmodels/sandbox/stats/multicomp.html

def ecdf(x):
    '''no frills empirical cdf used in fdrcorrection
    '''
    nobs = len(x)
    return np.arange(1,nobs+1)/float(nobs)

def fdrcorrection(pvals, alpha, method='indep'):

    '''pvalue correction for false discovery rate
    This covers Benjamini/Hochberg for independent or positively correlated and
    Benjamini/Yekutieli for general or negatively correlated tests. Both are
    available in the function multipletests, as method=`fdr_bh`, resp. `fdr_by`.
    Parameters
    ----------
    pvals : array_like
        set of p-values of the individual tests.
    alpha : float
        error rate
    method : {'indep', 'negcorr')
    Returns
    -------
    rejected : array, bool
        True if a hypothesis is rejected, False if not
    pvalue-corrected : array
        pvalues adjusted for multiple hypothesis testing to limit FDR
    Notes
    -----
    If there is prior information on the fraction of true hypothesis, then alpha
    should be set to alpha * m/m_0 where m is the number of tests,
    given by the p-values, and m_0 is an estimate of the true hypothesis.
    (see Benjamini, Krieger and Yekuteli)
    The two-step method of Benjamini, Krieger and Yekutiel that estimates the number
    of false hypotheses will be available (soon).
    Method names can be abbreviated to first letter, 'i' or 'p' for fdr_bh and 'n' for
    fdr_by.
    '''
    pvals = np.asarray(pvals)

    pvals_sortind = np.argsort(pvals)
    pvals_sorted = pvals[pvals_sortind]
    sortrevind = pvals_sortind.argsort()

    if method in ['i', 'indep', 'p', 'poscorr']:
        ecdffactor = ecdf(pvals_sorted)
    elif method in ['n', 'negcorr']:
        cm = np.sum(1./np.arange(1, len(pvals_sorted)+1))   #corrected this
        ecdffactor = ecdf(pvals_sorted) / cm
    elif method in ['n', 'negcorr']:
        cm = np.sum(np.arange(len(pvals)))
        ecdffactor = ecdf(pvals_sorted)/cm
    else:
        raise ValueError('only indep and necorr implemented')

    reject = pvals_sorted < ecdffactor*alpha
    if reject.any():
        rejectmax = max(np.nonzero(reject)[0])
    else:
        rejectmax = 0
    reject[:rejectmax] = True

    pvals_corrected_raw = pvals_sorted / ecdffactor
    pvals_corrected = np.minimum.accumulate(pvals_corrected_raw[::-1])[::-1]
    pvals_corrected[pvals_corrected>1] = 1
    return reject[sortrevind], pvals_corrected[sortrevind]
