#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import tensorflow as tf
from collections import UserDict


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


def split_train_test(dataset: tf.data.Dataset, test_prop: float = 0.3,
                     shuffle: bool = True, shuf_buf_size: int = 32,
                     batch_size: int = 16) -> (tf.data.Dataset, tf.data.Dataset):
    """Split a dataset into train and test dataset.

    Source: https://stackoverflow.com/questions/51125266/how-do-i-split-tensorflow-datasets
    """
    border = int(10 * test_prop)

    recover = lambda x, y: y
    is_test = lambda x, y: x % 10 < border
    is_train = lambda x, y: x % 10 >= border

    if shuffle:
        dataset = dataset.shuffle(shuf_buf_size)
    train_dtst = dataset.enumerate() \
            .filter(is_train) \
            .map(recover) \
            .batch(batch_size) \
            .cache()

    test_dtst = dataset.enumerate() \
            .filter(is_test) \
            .map(recover) \
            .batch(batch_size) \
            .cache()

    return train_dtst, test_dtst


def cmp(a, b):
    return (a > b) - (a < b)