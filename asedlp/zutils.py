#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
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
fmt = logging.Formatter("| {levelname: ^8} | {asctime} | {name}: {message}",
                        datefmt="%Y-%m-%d %H:%M:%S %p", style="{")
cs_handle = logging.StreamHandler()
cs_handle.setLevel(logging.DEBUG)
cs_handle.setFormatter(fmt)
logger.addHandler(cs_handle)
logger.setLevel(logging.DEBUG)


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


