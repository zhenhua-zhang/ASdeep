#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
from collections import UserDict

def print_arguments(arguments, fwidth=None, to_head=("subcommand", )):
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


class SmartDict(dict):
    """Container for the data set.
    """
    def __init__(self):
        super().__init__()
        
    def __setitem__(self, key, val):
        self.__dict__[key] = val

    def __getitem__(self, key):
        return self.__dict__[key]

    def __repr__(self):
        return repr(self.__dict__)

    def __delitem__(self, key):
        del self.__dict__[key]

    def clear(self):
        return self.__dict__.clear()

    def copy(self):
        return self.__dict__.copy()

    def aoc(self, key, val, create_list=True):
        """Append the value (if the key exists) or create the key (if not exists).
        """
        if key in self.__dict__:
            if isinstance(self.__dict__.get(key), list):
                self.__dict__[key].append(val)
            else:
                if create_list:
                    tmp_list = [self.__dict__[key]]
                    tmp_list.append(val)
                    self.__dict__[key] = tmp_list
                else:
                    raise TypeError("Require list, but found {}".format(type(self.__dict__.get(key))))
        else:
            self.__dict__[key] = val
