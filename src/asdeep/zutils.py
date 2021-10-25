# Author : Zhenhua Zhang
# E-mail : zhenhua.zhang217@gmail.com
# Created: May 07, 2020
# Updated: Oct 06, 2021

import re
import sys
import copy
import logging
import itertools as it

from collections import OrderedDict

import numpy as np
import matplotlib.pyplot as plt

# Reference: https://genomevolution.org/wiki/index.php/Ambiguous_nucleotide
M2B = { # monoallelic base-pair to biallelic base.
    "AA": "A", "AC": "M", "AG": "R", "AT": "W",
    "CA": "m", "CC": "C", "CG": "S", "CT": "Y",
    "GA": "r", "GC": "s", "GG": "G", "GT": "K",
    "TA": "w", "TC": "y", "TG": "k", "TT": "T",
    "NN": "N" }

B2M = { # Biallelic base to monoallelic base-pair
    "A": "AA", "M": "AC", "R": "AG", "W": "AT",
    "m": "CA", "C": "CC", "S": "CG", "Y": "CT",
    "r": "GA", "s": "GC", "G": "GG", "K": "GT",
    "w": "TA", "y": "TC", "k": "TG", "T": "TT",
    "N": "NN" }

# The length of chromosome by GRCh37. Ordered by Chromosome. Only autosomes.
# "X":  155270560, "Y":  59373566
CHRLEN = OrderedDict(
    {"1":  249250621, "2":  243199373, "3":  198022430, "4":  191154276,
     "5":  180915260, "6":  171115067, "7":  159138663, "8":  146364022,
     "9":  141213431, "10": 135534747, "11": 135006516, "12": 133851895,
     "13": 115169878, "14": 107349540, "15": 102531392, "16": 90354753,
     "17": 81195210,  "18": 78077248,  "19": 59128983,  "20": 63025520,
     "21": 48129895,  "22": 51304566})


# Logging manager.
class LogManager(logging.Logger):
    def __init__(self, name, level=logging.INFO, logstream: bool =True,
                 logfile: str=""):
        super(LogManager, self).__init__(name)

        fmt = logging.Formatter("{levelname: >8}|{asctime}|{name: >8}| {message}",
                                style="{", datefmt="%Y-%m-%d,%H:%M:%S")
        if logstream: self._add_hadler(logging.StreamHandler(), level, fmt)
        if logfile: self._add_hadler(logging.FileHandler(logfile), level, fmt)

    def _add_hadler(self, hdl, lvl, fmt):
        hdl.setLevel(lvl)
        hdl.setFormatter(fmt)
        self.addHandler(hdl)


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
            "`__getitem__()`, aka `[]`, methods and items()")

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
        img.save(fp=fp_out, format="GIF", append_images=imgs, save_all=True,
                duration=duration, loop=loop)
    except ValueError as err:
        print(err)


def fdr_bh(p):
    """A implementation of FDR by Benjamini & Hochberg.

    Reference:
        StackOverflow: https://stackoverflow.com/a/7451153/11203559
        BH = {
            i <- lp:1L   # lp is the number of p-values
            o <- order(p, decreasing = TRUE) # "o" will reverse sort the p-values
            ro <- order(o)
            pmin(1, cummin(n/i * p[o]))[ro]  # n is also the number of p-values }
    """
    p = np.asfarray(p)
    
    if (p < 0).any() or (1.0 < p).any():
        raise ValueError("P values should be between 0 and 1")

    l = len(p)
    o = np.argsort(-p)
    ro = np.argsort(o)
    
    return np.clip(np.minimum.accumulate(l / np.arange(l, 0, -1) * p[o])[ro], .0, 1.).tolist()


def fetch_layer_by_path(model, layer_path: str =".layer3[-1].conv1[0]"):
    regex = re.compile(r"\.\w+|\[-?\d+\]")
    per_access = re.findall(regex, layer_path)
    
    if len(per_access) == 0:
        return model

    for item in per_access:
        if item.startswith("."):
            model = getattr(model, item.replace(".", "", 1))
        elif item.startswith("["):
            model = model[int(item.replace("[", "").replace("]", ""))]
        else:
            raise ValueError("Unknown pattern was found: " + item)
    
    return model


def pickup_model(prebuilt_arch):
    """Pick up a modified models."""
    try:
        import torch.nn as nn
        import torchvision.models as models
    except ImportError as ime:
        logging.error("Import error: {}".format(ime))
        sys.exit(1)

    if prebuilt_arch == "resnext":
        # By default, using the pre-built ResNext.
        net = models.resnext50_32x4d(pretrained=False)
        net.conv1 = nn.Conv2d(1, 64, 11, 2, 3, bias=False)
        net.fc = nn.Linear(2048, 3, bias=True)
    elif prebuilt_arch == "resnet":
        net = models.resnet18(pretrained=False)
        net.conv1 = nn.Conv2d(1, 64, 11, 2, 3, bias=False)
        net.fc = nn.Linear(512, 3, bias=True)
    else:
        if prebuilt_arch != "alexnet":
            logging.warning("Unsupported prebuilt architecture, using default"
                            " AlexNet!")
        # Now I moved to AlexNet
        net = models.alexnet(pretrained=False)
        net.features[0] = nn.Conv2d(1, 64, 11, 4, 2)
        net.classifier[6] = nn.Linear(4096, 3, bias=True)

    return net


def parse_verbose(count):
    return min(60 - count * 10, 0)


def hcurve_to_img(hbcurve, output_prefix="./", cmap="summer", overlay=None,
                  overlay_cmap="Reds", overlay_alpha=0.5, scatter=True,
                  connect=True, color=True, kmer_text=True, figsize=(0, 0),
                  figfmt="pdf"):
    """Plot the Hilbert curve."""
    if sum(figsize) <= 0:
        _size = max(hbcurve.hcurve_matrix.shape) / 4
        figsize = (_size, _size)

    fig, ax = plt.subplots(figsize=figsize)

    if scatter:
        ax.scatter(hbcurve.x_pool, hbcurve.y_pool, c="0", marker=".", alpha=0.5)

    if connect:
        ax.plot(hbcurve.x_pool, hbcurve.y_pool, color="0", linewidth=0.5,
                alpha=0.5)

    if color:
        ax.imshow(hbcurve.hcurve_matrix, cmap=cmap)

        if isinstance(overlay, np.ndarray):
            ax.imshow(overlay, alpha=overlay_alpha, cmap=overlay_cmap)

    if kmer_text:
        for i, j in zip(hbcurve.x_pool, hbcurve.y_pool):
            i, j = int(i), int(j)
            _kmer_idx = int(hbcurve.hcurve_matrix[j, i])
            if _kmer_idx != -1:
                text = hbcurve.kmers_pool[_kmer_idx]
            else:
                text = "NULL"

            ax.text(i, j, text, ha="center", va="center",
                    fontsize="x-small", rotation=45)

    ax.set_axis_off()
    fig.set_tight_layout(True)

    save_path = "{}hilbert_curve.{}".format(output_prefix, figfmt)
    fig.savefig(save_path)

    plt.close(fig)


def calc_bits(seqlen, bits=0):
    if seqlen <= pow(2, bits)**2: return bits
    return calc_bits(seqlen, bits+1)

def make_all_mers(k, base="ACGT"):
    mers = it.product(base, repeat=k)
    mer2idx = {"N" * k: 0}
    mer2idx.update({"".join(m): i + 1 for i, m in enumerate(mers)})
    idx2mer = {v: k for k, v in mer2idx.items()}

    return mer2idx, idx2mer
