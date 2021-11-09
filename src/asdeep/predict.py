# Author : Zhenhua Zhang
# E-mail : zhenhua.zhang217@gmail.com
# Created: May 12, 2020
# Updated: Oct 06, 2021

"""Predict or validate the given input."""

import traceback
from argparse import Namespace
from collections import OrderedDict

import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.transforms as mtrans

import torch
from torch.nn import functional as func
from torch.autograd import Variable
from captum.attr import IntegratedGradients
from captum.attr import GradientShap
from captum.attr import Deconvolution
from captum.attr import GuidedGradCam

from .database import HilbertCurve
from .dataset import ASEDataset
from .dataset import SubsetHilbertCurve
from .zutils import fetch_layer_by_path
from .zutils import pickup_model
from .zutils import LogManager


class MissingSampleIdError(Exception):
    pass


class NegativeLabelError(Exception):
    pass


class NonSupportedMatrixType(Exception):
    pass


class Predictor:
    """Predict the new sample using a given model."""
    def __init__(self, model, model_state, dataset: ASEDataset,
                 model_arch: str = "alexnet", store_attrs: bool = False,
                 logman: LogManager = LogManager("Predict")):
        self._logman = logman
        # Whether store attributions
        self._store_attrs = store_attrs

        # Check if CUDA is available
        self._device = "cuda:0" if torch.cuda.is_available() else "cpu"

        # Load the model from given model state file or a model
        self._model = self._load_model(model, model_state, model_arch)

        # The database contains Hilbert curve converted from sequence
        self._dataset = dataset

        # Containers to store prediction results and sample information
        self._results = OrderedDict()

        # Container to store attributions
        self._attrs = {}

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, tb):
        if exc_type is not None:
            traceback.print_exception(exc_type, exc_value, tb)

        self

    def _load_model(self, model, model_state, model_arch):
        state_dict = torch.load(model_state, map_location=self.device)
        model.load_state_dict(state_dict[model_arch])
        model.eval()
        return model

    def _load_sample(self, sample_id):
        if sample_id in self._dataset:
            hbcurve, attrs = self._dataset[sample_id]
            return (hbcurve, attrs)
        raise MissingSampleIdError("Sample ID was not found!")

    def _predict(self, matrix):
        # Predict the given dataset.
        if isinstance(matrix, HilbertCurve):
            matrix = matrix.hbcmat

        matrix = Variable(torch.Tensor(np.expand_dims(matrix, 0)))
        logit = self._model(matrix)
        prob, label = func.softmax(logit, dim=1).data.squeeze().sort(0, True)
        return prob.numpy()[0], label.numpy()[0]

    def _calc_attrs(self, matrix, label: int, attr_mthd: str,
                    gdc_layer_path: str=".layer3[-1].conv3"):
        """Calculate attributions."""
        # Fetch matrix
        if isinstance(matrix, HilbertCurve):
            matrix = matrix.hbcmat

        matrix = Variable(torch.Tensor(np.expand_dims(matrix, 1)))

        label = int(label) # Force label as int.
        if label < 0:
            raise NegativeLabelError("Labels should be zero or non-negative")

        if attr_mthd in ("IG", "IntegratedGradients"): # IntegratedGradients
            attr = (IntegratedGradients(self._model)
                    .attribute(matrix, target=label, n_steps=200))
        elif attr_mthd in ("GS", "GradientShap"): # GradientShap
            rand_mtx_dist = torch.cat([matrix*0, matrix*1, matrix*2])
            attr = (GradientShap(self._model)
                    .attribute(matrix, target=label, n_samples=50,
                                  baselines=rand_mtx_dist))
        elif attr_mthd in ("GGC", "GuidedGradCam"): # GuidedGradCam
            gradcam_layer = fetch_layer_by_path(self._model, gdc_layer_path)
            attr = (GuidedGradCam(self._model, gradcam_layer)
                    .attribute(matrix, target=label))
        else: # Deconvolution, the default
        # NOTE: Potential issue/warnings
        #     1. https://github.com/pytorch/captum/issues/408
            attr = Deconvolution(self._model).attribute(matrix, target=label)
            if attr_mthd not in ("DC", "Deconvolution"):
                self._logman.warning(f"Unsupported attribute: {attr_mthd}."
                                     " Using Deconvolution by default.")

        return np.expand_dims(attr.squeeze().cpu().detach().numpy(), 0)

    @property
    def model(self):
        return self._model

    @property
    def device(self):
        return self._device

    @property
    def predictions(self):
        return {s: {"sample": s, "label": l, "prob": p}
                for s, (l, p, _) in self._results.items()}

    @property
    def attributions(self):
        return self._attrs

    def predict(self, sample_ids: list, keep_attrs=None):
        """Predict the given input.

        Args:
        Returns:
        Raise:
        """
        # _label = {0: "ASEto0", 1: "NonASE", 2: "ASEto2", None: "Unknown"}
        self._logman.info("Sample ID, Pred label, Pred prob")
        for per_sample in sample_ids:
            hbcurve, hbcattr = self._load_sample(per_sample)
            prob, label = self._predict(hbcurve)
            self._logman.info(f"{per_sample: >9}, {label: >10}, {prob: >10.3}")

            self._results.update({per_sample: [label, prob, hbcattr]})

            if keep_attrs is None or len(keep_attrs) < 1:
                continue

            if "all" in keep_attrs:
                keep_attrs = ["IG", "GS", "DC"]#, "GGC"]

            attr_p_sample = OrderedDict({"HBC": hbcurve})
            for p_attr in keep_attrs:
                if p_attr not in self._attrs:
                    p_attr_hbc = self._calc_attrs(hbcurve, label, p_attr)
                    attr_p_sample[p_attr] = p_attr_hbc
                else:
                    self._logman.warning(f"{p_attr} alreay added, skip it.")

            self._attrs[per_sample] = attr_p_sample

        return self

    def show_attrs(self, figsize=(16, 9), scale=1000, output_dir="./",
                   fmt="png"):
        """Show the attributions.

        Args:
        Returns:
        Raise:
        """
        for p_sample, p_attr_map in self._attrs.items():
            nattrs = len(p_attr_map)

            sample_info = self._results[p_sample]
            strand = sample_info[-1].get("strand", 1)

            fig_width, fig_height = figsize
            mat_ratio = int(fig_height/nattrs)
            seq_ratio = fig_width - mat_ratio
            gridspec_kw={"width_ratios": [mat_ratio, seq_ratio]}

            if nattrs < 1:
                self._logman.warning("No attribution were calculated.")
            else:
                fig, axe = plt.subplots(nrows=nattrs, ncols=2, figsize=figsize,
                                        gridspec_kw=gridspec_kw)

                self._logman.info(f"Plotting for {p_sample}")
                for idx, (p_attr, p_attr_mat) in enumerate(p_attr_map.items()):
                    if isinstance(p_attr_mat, HilbertCurve):
                        attr_hbc = p_attr_mat
                    else:
                        attr_hbc = HilbertCurve(p_attr_mat, strand=strand,
                                                dtype=np.float64)

                    attr_hbc_mat = attr_hbc.get_hbcmat(0.0).squeeze() * scale
                    attr_a1_seq, attr_a2_seq = attr_hbc.allelic_attrs

                    norm = mcolors.TwoSlopeNorm(vcenter=0)
                    base = axe[idx, 0].transData
                    tran = mtrans.Affine2D().rotate_deg(90)
                    axe[idx, 0].pcolormesh(attr_hbc_mat, cmap="seismic",
                                           norm=norm, transform=tran+base)

                    axe[idx, 0].set_xticks([])
                    axe[idx, 0].set_yticks([])

                    attr_x_vals = np.arange(len(attr_a1_seq))
                    if p_attr == "HBC":
                        is_het = attr_a2_seq != attr_a1_seq
                        het_pos = attr_x_vals[is_het]
                        het_site = (np.zeros_like(is_het) + 1)[is_het]
                        axe[idx, 1].scatter(het_pos, het_site)
                        axe[idx, 1].set_yticks([1])
                        axe[idx, 1].set_ylabel(f"Heterozygous SNPs")
                    else:
                        attr_allelic_diff = attr_a1_seq - attr_a2_seq
                        ylim = max(abs(attr_allelic_diff))
                        axe[idx, 1].plot(attr_x_vals, attr_allelic_diff)
                        axe[idx, 1].set_ylabel(f"Attribution diff ({p_attr})")
                        axe[idx, 1].set_ylim((-ylim, ylim))

                    if idx < nattrs - 1:
                        axe[idx, 1].spines["bottom"].set_visible(False)
                        axe[idx, 1].set_xticks([])

                    axe[idx, 1].spines["top"].set_visible(False)
                    axe[idx, 1].spines["left"].set_visible(False)
                    axe[idx, 1].set_xlim((0, len(attr_x_vals)))
                    axe[idx, 1].yaxis.tick_right()

                fig.set_tight_layout(True)
                fig.savefig(f"{output_dir}/{p_sample}_attributions.{fmt}")

        return self


def predict(args: Namespace, logman: LogManager=LogManager("Predict")):
    """Predict."""
    database = args.database
    sample_ids = args.sample_ids
    output_dir = args.output_dir
    attributions = args.attributions
    prebuilt_arch = args.prebuilt_arch
    model_state_path = args.model_state_path
    save_fmt = args.save_fmt

    logman.info(f"Database:     {database}")
    logman.info(f"Sample IDs:   {sample_ids}")
    logman.info(f"Model path:   {model_state_path}")
    logman.info(f"Attribution:  {attributions}")
    logman.info(f"Output path:  {output_dir}")
    logman.info(f"Architecture: {prebuilt_arch}")

    net = pickup_model(prebuilt_arch)
    trans = [SubsetHilbertCurve()]
    with h5.File(database, "r") as h5db:
        dataset = ASEDataset(database=h5db, transformers=trans)
        with Predictor(net, model_state_path, dataset, prebuilt_arch) as pred:
            (pred
             .predict(sample_ids, attributions)
             .show_attrs(output_dir=output_dir, fmt=save_fmt))
