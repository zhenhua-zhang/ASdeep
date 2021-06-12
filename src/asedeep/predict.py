#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Predict or validate the given input.
"""

import logging
from typing import Union, List, Tuple, Iterable

import numpy as np
import matplotlib.pyplot as plt

import torch
import torch.nn as nn
from torch.nn import functional as func
from torch.autograd import Variable
from torchvision import models
from captum.attr import IntegratedGradients
# from captum.attr import GradientShap
# from captum.attr import Deconvolution
# from captum.attr import GuidedGradCam

from .dataset import ASEDataset
from .dataset import SeqToHilbertAndMakeLabel
from .dataset import MultipleTestAdjustment
from .hbcurve import HilbertCurve
from .zutils import fetch_layer_by_path


class Predictor(object):

    def __init__(self,
                 net: torch.nn.Module=None,
                 gene_id: str=None,
                 file_path: str=None):
        self.net = net
        self.gene_id = gene_id
        self.file_path = file_path

        self._dataset: ASEDataset

        self._matrix: np.ndarray
        self._hbcurve: HilbertCurve
        self._sequence: np.ndarray

        self._pred_prob: float
        self._pred_label: int
        self._true_label: int

        self._attrseq: np.ndarray
        self._attrmtx: np.ndarray

    @staticmethod
    def _check_device():
        return "cuda:0" if torch.cuda.is_available() else "cpu"

    def _check_matrix_idx(self,
                          matrix_idx):
        if matrix_idx is None:
            return range(len(self._dataset))
        elif isinstance(matrix_idx, int):
            return [matrix_idx]
        elif isinstance(matrix_idx, (list, tuple)) \
                or hasattr(matrix_idx, "__next__"):
            return matrix_idx
        else:
            raise TypeError("Unsupported matrix_idx: " + type(matrix_idx))

    def _load_mtx(self,
                  matrix_idx: Union[int, List[int], Tuple[int], Iterable[int]]):
        self._hbcurve, (_, self._true_label) = self._dataset[matrix_idx]
        self._matrix = self._hbcurve.get_hcurve()
        self._sequence = self._hbcurve.get_kmerseq()

    def _predict(self):
        # Predict the given dataset.
        matrix = Variable(torch.Tensor(np.expand_dims(self._matrix, 0)))
        logit = self.net(matrix)
        probs, idx = func.softmax(logit, dim=1).data.squeeze().sort(0, True)
        self._pred_prob, self._pred_label = probs.numpy()[0], idx.numpy()[0]

    def _attrseq_to_nucseq(self):
        """Convert attribution sequence into nucleotide sequence."""
        self._dataset

    def _calc_attrmtx(self,
                      gdc_layer_path: str=".layer3[-1].conv3"):
        """Calculate attributions.
        """
        matrix = Variable(torch.Tensor(np.expand_dims(self._matrix, 0)))
        pred_label = int(self._pred_label)
        # IntegratedGradients
        attr_ig = IntegratedGradients(self.net).attribute(matrix, target=pred_label,
                                                          n_steps=200)

        # # GradientShap
        # rand_mtx_dist = torch.cat([matrix*0, matrix*1, matrix*2])
        # attr_gs = GradientShap(self.net).attribute(matrix, target=pred_label, n_samples=50,
        #                                            baselines=rand_mtx_dist)

        # # Occlusion
        # attr_oc = Occlusion(self.net).attribute(matrix, target=pred_label, baselines=0,
        #                                         strides=(1, 2, 2), sliding_window_shapes=(1, 4, 4))

        # # Deconvolution
        # attr_dc = Deconvolution(self.net).attribute(matrix, target=pred_label)

        # # GuidedGradCam
        # gradcam_layer = fetch_layer_by_path(self.net, gdc_layer_path)
        # attr_gd = GuidedGradCam(self.net, gradcam_layer).attribute(matrix, target=pred_label)

        attr_method = ["Integrated gradients"] # ["Gradient SHAP", "Deconvolution", "Guided gradient CAM"]
        attrs = [attr_ig] # [attr_gs, attr_dc, attr_gd]

        self._attrmtx = np.zeros(self._matrix.shape[-2:], dtype=[(k, ">f8") for k in attr_method])
        for mthd, attr in zip(attr_method, attrs):
            self._attrmtx[mthd] = attr.squeeze().cpu().detach().numpy()

    def _plot_attrmtx(self,
                      orig_color="Greys",
                      attr_color="Reds",
                      figtitle="Hilbert Curve",
                      figsize=(8, 8)):
        """Draw heatmaps to show the attribution of a few approaches.

        Args:
            orig_color:
            attr_color:
            figtitle:
            figsize:

        Returns:
            A tuple of (Figure, Axes), aka a Figure instance and an Axes instance.
        """
        fig, axes = plt.subplots(2, len(self._attrmtx.dtype) + 1, figsize=figsize)

        axes[0][0].annotate("Attribution per DNA kmer.", (.5, .5), ha="center",
                            fontsize="x-large")
        axes[0][0].set_axis_off()

        axes[1][0].imshow(self._matrix.squeeze(), alpha=1.0, cmap=orig_color)
        axes[1][0].set_axis_off()
        axes[1][0].set_title(figtitle)

        for idx, mthd in enumerate(self._attrmtx.dtype.names):
            ax_idx = idx + 1
            axes[0][ax_idx].imshow(self._attrmtx[mthd], alpha=1.0, cmap=attr_color)
            axes[0][ax_idx].set_axis_off()
            axes[0][ax_idx].set_title(mthd)

            axes[1][ax_idx].imshow(self._matrix.squeeze(), alpha=1.0, cmap=orig_color)
            axes[1][ax_idx].imshow(self._attrmtx[mthd], alpha=0.5, cmap=attr_color)
            axes[1][ax_idx].set_axis_off()

        fig.set_tight_layout(tight=True)
        
        return (fig, axes)

    def _attrmtx_to_attrseq(self,
                            scale_hbcv: bool =True,
                            biallelic: bool =True):
        """Converts attributions of the Hilbert curve matrix into a sequence.
        
        Args:
            scale_hbcv:
            biallelic:

        Returns:
            A named numpy.narray of attributions from given attribution matrix. The
            results include all attributions given by attr_map as well as the code
            of Hilbert curve matrix.

        Raises:
            None
        """
        if scale_hbcv: self._matrix *= 1.0 / self._matrix.max()
        
        n_elements = np.prod(self._matrix.shape)  # Number of elements in the matrix
        _attrseq = np.zeros(n_elements, dtype=self._attrmtx.dtype)
        for idx, (i, j) in enumerate(zip(self._hbcurve.x_pool, self._hbcurve.y_pool)):
            i, j = int(i), int(j)
            for key in self._attrmtx.dtype.names:
                _attrseq[key][idx] = self._attrmtx[key][j, i]

        if biallelic:
            allele_len = int(n_elements / 2)

            # The original sequence
            _sequence = self._sequence
            self._sequence = np.zeros((allele_len, 2))
            self._sequence[:, 0] = _sequence[:allele_len]
            self._sequence[:, 1] = _sequence[allele_len:][::-1]

            # The attributions
            self._attrseq = np.zeros((allele_len, 2), dtype=self._attrmtx.dtype)
            for key in _attrseq.dtype.names:
                self._attrseq[key][:, 0] = _attrseq[key][:allele_len]
                self._attrseq[key][:, 1] = _attrseq[key][allele_len:][::-1]
        else:
            self._attrseq = _attrseq
            
    def _plot_attrseq(self,
                      colors: tuple =("r", "b"),
                      figsize: tuple =(16, 8)):
        """Draw a line plot to show the attribution of the n-mer along the sequence.

        Args:
            attrs:
                A numpy.ndarray of attributions with type of attribution as name and 
                attributions sequence as value. The attribution sequences should be
                splited into two alleles.
            colors:
                A tuple of colors including two matplotlib.color, one for allele A
                the other for allele B
            figsize:
                A tuple of floats to set the figure size.

        Returns:
            A tuple including matplotlib.pyplot.Figure and matplotlib.pyplot.axes.

        Raises:
            ValueError:
                A error occur when the sequence of Hilbert curve not splited into
                sequences of two allele.
        """
        seqlen, n_allele = self._attrseq.shape
        if n_allele != 2:
            raise ValueError("The sequence should be splited into two allele!")

        # Check if there is any heterogeneous sites.
        heter_sites = (self._sequence[:, 0] != self._sequence[:, 1]).nonzero()
        if heter_sites[0].size < 1:
            heter_sites = None
        else:
            heter_sites = heter_sites[0]            

        # plot sequence
        attr_names = self._attrseq.dtype.names

        n_attrs = len(attr_names)
        fig, axes = plt.subplots((n_attrs + 1) * 2 + 1, figsize=figsize, sharex=True)

        x_coord = np.arange(seqlen)
        for idx, attr_name in enumerate(attr_names):
            for axes_idx in [idx * 2, idx * 2 + 1]:
                attr_idx = axes_idx - idx * 2
                axes[axes_idx].plot(x_coord, self._attrseq[attr_name][:, attr_idx],
                                    color=colors[attr_idx])
                axes[axes_idx].set_title("{} (allele {})".format(attr_name, attr_idx))

                axes[axes_idx].set_xticks([])        
                axes[axes_idx].spines["top"].set_visible(False)
                axes[axes_idx].spines["right"].set_visible(False)
                axes[axes_idx].spines["bottom"].set_visible(False)

        for axes_idx in [-2, -3]:
            seq_idx = axes_idx + 3
            axes[axes_idx].plot(x_coord, self._sequence[:, seq_idx], color=colors[seq_idx])
            axes[axes_idx].set_title("Original sequence (allele {})".format(seq_idx))
            axes[axes_idx].set_xticks([])        
            axes[axes_idx].spines["top"].set_visible(False)
            axes[axes_idx].spines["right"].set_visible(False)
            axes[axes_idx].spines["bottom"].set_visible(False)

        if heter_sites is not None and heter_sites.size > 0:
            axes[-1].vlines(heter_sites, 0, 1, color="k")
        axes[-1].set_xticks([0, max(x_coord)])
        axes[-1].set_title("Heterogeneous sites")
        axes[-1].spines["top"].set_visible(False)
        axes[-1].spines["right"].set_visible(False)

        fig.set_tight_layout(True)
        fig.align_ylabels()

        return (fig, axes)

    def init(self,
             model_state,
             output_size=64):
        """Initialize to have basic settings ready.

        Args:
            model_state:
            output_size:

        Returns:
            The Predictor it self.

        Raises:
            None
        """
        device = self._check_device()

        if self.net is None: # Load model
            self.net = models.resnext50_32x4d()

            # in_channels 1, out_channels 64, kernel_size 7, stride 2, padding 3
            self.net.conv1 = nn.Conv2d(1, output_size, 7, 2, 3, bias=False)
            self.net.fc = nn.Linear(2048, 3)

        self.net.load_state_dict(torch.load(model_state, map_location=device))
        self.net.eval()

        return self

    def load_dataset(self,
                     gene_id=None,
                     file_path_pool=None,
                     element_trans=None,
                     dataset_trans=None):
        """Load dataset.

        Args:
            gene_id:
            file_path_pool:
            element_trans:
            dataset_trans:

        Returns:
            The instance of Predictor itself.
        """
        if gene_id is None:
            gene_id = self.gene_id

        if file_path_pool is None:
            file_path_pool = self.file_path

        if element_trans is None:
            element_trans = SeqToHilbertAndMakeLabel(matrix=False)

        if dataset_trans is None:
            dataset_trans = MultipleTestAdjustment()

        self._dataset = ASEDataset(file_path_pool=file_path_pool, gene_id=gene_id,
                                   element_trans=element_trans, dataset_trans=dataset_trans)

        return self

    def predict(self,
                save_pref,
                matrix_idx=None,
                show_attr=False):
        """Predict the given input.

        Args:
            save_pref:
            matrix_idx:
            show_attr:

        Returns:
            The instance of Predictor itself.
        """
        matrix_idx = self._check_matrix_idx(matrix_idx)
        label = {0: "ASEtoA", 1: "NonASE", 2: "ASEtoB", None: "Unknown"}

        logging.info("TrueLabel | PredLabel | PredProb")
        for idx in matrix_idx:
            self._load_mtx(idx)
            self._predict()
            true_label, pred_label = label[self._true_label], label[self._pred_label]
            logging.info("{: ^9} | {: ^9} | {: ^9.3}".format(true_label, pred_label,
                                                             self._pred_prob))

            if show_attr:
                splitter = "" if save_pref.endswith("/") else "."
                _save_name = "{}{}{}_{}_{}.pdf".format(save_pref, splitter, idx, true_label,
                                                       pred_label)
                self._calc_attrmtx()
                fig, _ = self._plot_attrmtx()
                fig.savefig(_save_name.replace(".pdf", "_attrmtx.pdf"))

                self._attrmtx_to_attrseq()
                fig, _ = self._plot_attrseq()
                fig.savefig(_save_name.replace(".pdf", "_attrseq.pdf"))

        return self


if __name__ == "__main__":
    logging.warning("This module should not be executed directly.")
