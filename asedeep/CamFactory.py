#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""An implementation of Class Activation Mapping (CAM) based on 

Reference: https://github.com/zhoubolei/CAM/blob/master/pytorch_CAM.py
"""

import cv2
import torch
import numpy as np
import torch.nn as nn
import matplotlib.pyplot as plt

from torchvision import models
from torch.autograd import Variable
from torch.nn import functional as func

from .ASEDataset import ASEDataset, ReshapeMatrixAndPickupLabel, MultipleTestAdjustMent


class CamFactory:
    def __init__(self, net=None, gene_id=None, file_path=None, feature_layer="layer4"):
        self.net = net
        self.gene_id = gene_id
        self.file_path = file_path
        self.feature_layer = feature_layer

        self.feature_blobs = []
        self.label = None
        self.matrix = None
        self.pred_prob = None
        self.pred_label = None

        self.mtx_img = None
        self.dataset = None
        self.cam_img = None
        self.cam_hmap = None
        self.seq_img = None
        self.seq_hmap = None
        self.weight_softmax = None

    @staticmethod
    def _check_device():
        return "cuda:0" if torch.cuda.is_available() else "cpu"

    def _check_matrix_idx(self, matrix_idx):
        if matrix_idx is None:
            return range(len(self.dataset))
        elif isinstance(matrix_idx, int):
            return [matrix_idx]
        elif isinstance(matrix_idx, (list, tuple)) or hasattr(matrix_idx, "__next__"):
            return matrix_idx
        else:
            raise TypeError("Unsupported matrix_idx: {}.".format(type(matrix_idx)))

    def _reg_hook(self, hook=None, flayer="layer4"):
        if flayer not in self.net._modules:
            raise KeyError("Uknown feature layer: {}".format(flayer))

        def feature_hook(module, input, output):
            self.feature_blobs.append(output.data.cpu().numpy())

        if hook is None:
            hook = feature_hook

        self.net._modules.get(flayer).register_forward_hook(hook)

    def _load_matrix(self, matrix_idx):
        matrix, self.label = self.dataset[matrix_idx]
        self.matrix = np.expand_dims(matrix, axis=0)

    def _mtx_to_img(self, matrix=None):
        if matrix is None:
            matrix = self.matrix

        matrix = matrix.squeeze(0)
        c, w, h = matrix.shape
        self.mtx_img = np.uint8(matrix.reshape((w, h, c)) * 255)
        # pdb.set_trace()

    def _calc_cam(self):
        # Calculate class activation map.
        matrix = Variable(torch.Tensor(self.matrix))
        logit = self.net(matrix)
        prob_and_idx = func.softmax(logit, dim=1).data.squeeze()
        probs, idx = prob_and_idx.sort(0, True)
        probs, idx = probs.numpy(), idx.numpy()
        self.pred_prob, self.pred_label = probs[0], idx[0]
        # pdb.set_trace()

        feature_conv = self.feature_blobs[-1]
        bz, nc, h, w = feature_conv.shape
        cam = self.weight_softmax[idx[0]].dot(feature_conv.reshape((nc, h*w)))
        cam = cam.reshape(h, w)
        cam = cam - np.min(cam)
        cam = cam / np.max(cam)
        cam_img = np.uint8(255 * cam)

        self.cam_img = cv2.resize(cam_img, (self.matrix.shape[-2:]))
        self.cam_hmap = cv2.applyColorMap(self.cam_img, cv2.COLORMAP_JET)

    def _cam_to_seq(self):
        width, height, channel = self.cam_hmap.shape
        self.seq_hmap = self.cam_hmap.reshape(int(width * height / 8), 8, channel)
        # pdb.set_trace()

        width, height = self.cam_img.shape
        self.seq_img = self.cam_img.reshape(int(width * height / 8), 8)

    def _onehot_to_text(self):
        # Placeholder to implement function convert one-hot encoded sequence
        # into nucleotide sequence.
        self.seq_img

    def init(self, model_state):
        device = self._check_device()

        # Load model
        if self.net is None:
            self.net = models.resnext50_32x4d()
            # in_channels 1, out_channels 64, kernel_size 7, stride 2, padding 3
            self.net.conv1 = nn.Conv2d(1, 64, 7, 2, 3, bias=False)
            self.net.fc = nn.Linear(2048, 3)

        self.net.load_state_dict(torch.load(model_state, map_location=device))
        self.net.eval()

        # Fetch parameters
        params = list(self.net.parameters())
        self.weight_softmax = np.squeeze(params[-2].data.numpy())

        self._reg_hook(flayer=self.feature_layer)

        return self

    def load_dataset(self, **kwargs):
        """Load dataset.
        """
        gene_id = kwargs.get("gene_id", self.gene_id)
        file_path_pool = kwargs.get("file_path_pool", self.file_path)
        element_trans = kwargs.get("element_trans", ReshapeMatrixAndPickupLabel())
        dataset_trans = kwargs.get("dataset_trans", MultipleTestAdjustMent())

        self.dataset = ASEDataset(file_path_pool=file_path_pool, gene_id=gene_id, element_trans=element_trans, dataset_trans=dataset_trans)

        return self

    def show_cam(self, save_path, matrix_idx=None):
        matrix_idx = self._check_matrix_idx(matrix_idx)
        for idx in matrix_idx:
            self._load_matrix(idx)
            self._mtx_to_img()
            self._calc_cam()
            self._cam_to_seq()

            cam_hmap = self.mtx_img * 0.5 + self.cam_hmap * 0.5
            cv2.imwrite(save_path + ".{}_cam_max_hmap_{}.png".format(idx, self.label), cam_hmap)
            cv2.imwrite(save_path + ".{}_mtx_hmap_{}.png".format(idx, self.label), self.mtx_img)
            cv2.imwrite(save_path + ".{}_cam_hmap_{}.png".format(idx, self.label), self.cam_hmap)

        return self

    def _draw_cam_along_seq(self, save_path, title):
        fig, ax = plt.subplots()
        height, width = self.seq_img.shape
        average_cam = np.sum(self.seq_img, axis=1) / width / 255    # 255 is the scale factor

        ax.plot(range(height)[::-1], average_cam)
        ax.set_title(title)
        fig.savefig(save_path)
        plt.close(fig)

    def show_cam_dist_along_seq(self, save_path, matrix_idx=None):
        matrix_idx = self._check_matrix_idx(matrix_idx)

        lable_dict = {0: "ASEtoA", 1: "NonASE", 2: "ASEtoB"}

        for idx in matrix_idx:
            self._load_matrix(idx)
            self._mtx_to_img()
            self._calc_cam()
            self._cam_to_seq()

            pic_save_name = ".{}.cam_along_seq.{}.{}.png".format(idx, self.label, self.pred_label)
            true_label, pred_label = lable_dict[self.label], lable_dict[self.pred_label]
            pic_title = "sample-{:0>3}: {}(true label), {}(pred label)".format(idx, true_label, pred_label)
            self._draw_cam_along_seq(save_path + pic_save_name, pic_title)

        return self
