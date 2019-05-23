#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys

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
