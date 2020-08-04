#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# File name : zutils.py
# Author    : Zhenhua Zhang
# E-mail    : zhenhua.zhang217@gmail.com
# Created   : Mon 22 Jun 2020 11:25:03 AM CEST
# Version   : v0.1.0
# License   : MIT
#

import sys
import copy
import logging

from collections import UserDict

logger = logging.getLogger()
fmt_str = "| {levelname: ^8} | {asctime} | {name}: {message}"
dt_str = "%Y-%m-%d %H:%M:%S %p"
fmt = logging.Formatter(fmt_str, datefmt=dt_str, style="{")
cs_handle = logging.StreamHandler()
cs_handle.setLevel(logging.INFO)
cs_handle.setFormatter(fmt)
logger.addHandler(cs_handle)
logger.setLevel(logging.INFO)


# Reference: https://genomevolution.org/wiki/index.php/Ambiguous_nucleotide
# M2B = { # monoallelic base-pair to biallelic base.
#     "AA": "A", "AC": "M", "AG": "R", "AT": "W",
#     "CA": "m", "CC": "C", "CG": "S", "CT": "Y",
#     "GA": "r", "GC": "s", "GG": "G", "GT": "K",
#     "TA": "w", "TC": "y", "TG": "k", "TT": "T",
#     "NN": "N"
# }

M2B = { # monoallelic base-pair to biallelic base.
    (65, 65): "A", (65, 67): "M", (65, 71): "R", (65, 84): "W",
    (67, 65): "m", (67, 67): "C", (67, 71): "S", (67, 84): "Y",
    (71, 65): "r", (71, 67): "s", (71, 71): "G", (71, 84): "K",
    (84, 65): "w", (84, 67): "y", (84, 71): "k", (84, 84): "T",
    (78, 78): "N",
}

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
#                      shuffle: bool = True, shuf_buf_size: int = 32,
#                      batch_size: int = 16) -> (tf.data.Dataset, tf.data.Dataset):
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



