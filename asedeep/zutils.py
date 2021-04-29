#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author     : Zhenhua Zhang
# E-mail     : zhenhua.zhang217@gmail.com
# Created on : Mon 22 Jun 2020 11:25:03 AM CEST


import sys
import copy
import numpy as np

from collections import UserDict, OrderedDict

# Reference: https://genomevolution.org/wiki/index.php/Ambiguous_nucleotide
M2B = { # monoallelic base-pair to biallelic base.
    'AA': 'A', 'AC': 'M', 'AG': 'R', 'AT': 'W',
    'CA': 'm', 'CC': 'C', 'CG': 'S', 'CT': 'Y',
    'GA': 'r', 'GC': 's', 'GG': 'G', 'GT': 'K',
    'TA': 'w', 'TC': 'y', 'TG': 'k', 'TT': 'T',
    'NN': 'N'
}

B2M = { # Biallelic base to monoallelic base-pair
    'A': 'AA', 'M': 'AC', 'R': 'AG', 'W': 'AT',
    'm': 'CA', 'C': 'CC', 'S': 'CG', 'Y': 'CT',
    'r': 'GA', 's': 'GC', 'G': 'GG', 'K': 'GT',
    'w': 'TA', 'y': 'TC', 'k': 'TG', 'T': 'TT',
    'N': 'NN'
}

# The length of chromosome by GRCh37. Ordered by Chromosome. Only autosomes.
# 'X':  155270560, 'Y':  59373566
CHRLEN = OrderedDict(
    {'1':  249250621, '2':  243199373, '3':  198022430, '4':  191154276, '5':  180915260,
     '6':  171115067, '7':  159138663, '8':  146364022, '9':  141213431, '10': 135534747,
     '11': 135006516, '12': 133851895, '13': 115169878, '14': 107349540, '15': 102531392,
     '16': 90354753,  '17': 81195210,  '18': 78077248,  '19': 59128983,  '20': 63025520,
     '21': 48129895,  '22': 51304566})


def print_arguments(arguments, fwidth: int=0, to_head=('subcommand',)):
    '''Print artuments from command lines'''
    print('Arguments for current run: ', file=sys.stderr)
    _dict = vars(arguments)
    _pair = [(_dst, _arg) for _dst, _arg in _dict.items() if _dst not in to_head]
    _pair = sorted(_pair, key=lambda x: len(x[0]))

    _to_head = [(_dst, _arg) for _dst, _arg in _dict.items() if _dst in to_head]

    for prior in _to_head:
        _pair.insert(0, prior)

    if fwidth is None:
        fwidth = len(_pair[-1][0]) + 1

    for _dst, _arg in _pair:
        if _arg is None:
            _arg = 'None'

        if isinstance(_arg, (list, tuple)):
            _arg = ', '.join(map(str, _arg))

        print('    {d: <{w}}: {a: <{w}}'.format(d=_dst, a=_arg, w=fwidth), file=sys.stderr)


class SmartDict(UserDict):
    '''Container for the data set.
    '''

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
    '''Flatten a list recursively.
    Args:
        l (list): The list to be flatten.
    Returns:
        out_list (list): The flatten list.
    '''
    out_list = []
    for x in l:
        if isinstance(x, (list, tuple)):
            out_list.extend(flatten(x))
        else:
            out_list.append(x)
    return out_list


def insert_or_append(dict1, dictlike2):
    dict1 = copy.deepcopy(dict1)

    if not hasattr(dictlike2, '__getitem__') or not hasattr(dictlike2, 'items'):
        raise AttributeError(
            'dictlike2 should at least has '
            '`__getitem__()`, aka `[]`, methods and items()')

    for key, val in dictlike2.items():
        if key in dict1:
            dict1[key].append(val)
        else:
            dict1[key] = [val]

    return dict1


def make_gif(fp_in, fp_out, duration=400, loop=0):
    import glob
    from PIL import Image

    try:
        fp_in_imgs = glob.glob(fp_in)
        img, *imgs = [Image.open(f) for f in sorted(fp_in_imgs)]
        img.save(fp=fp_out, format='GIF', append_images=imgs, save_all=True,
                duration=duration, loop=loop)
    except ValueError as err:
        print(err)


def fdr_bh(p):
    '''A implementation of FDR by Benjamini & Hochberg.

    Reference:
        StackOverflow: https://stackoverflow.com/a/7451153/11203559
        BH = {
            i <- lp:1L   # lp is the number of p-values
            o <- order(p, decreasing = TRUE) # "o" will reverse sort the p-values
            ro <- order(o)
            pmin(1, cummin(n/i * p[o]))[ro]  # n is also the number of p-values }
    '''
    p = np.asfarray(p)
    
    if (p < 0).any() or (1.0 < p).any():
        raise ValueError("P values should be between 0 and 1")

    l = len(p)
    o = np.argsort(-p)
    ro = np.argsort(o)
    
    return np.clip(np.minimum.accumulate(l / np.arange(l, 0, -1) * p[o])[ro], .0, 1.).tolist()

