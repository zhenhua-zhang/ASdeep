#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""An implementation of Class Activation Mapping (CAM) based on 

Reference: https://github.com/zhoubolei/CAM/blob/master/pytorch_CAM.py
"""

import logging
import sys

import cv2
import numpy as np

import torch
from torch.nn import functional as func

from .zutils import logger


class CamFactory:
    def __init__(self, net, matrix, feature_layer="layer4"):
        self.net = net
        self.matrix = matrix
        self.feature_layer = feature_layer

        self.feature_blobs = []
        self.cam_img = None
        self.heatmap = None

    def _reg_hook(self, feature_layer="layer4"):
        if feature_layer not in self.net._modules:
            raise KeyError("Uknown feature layer: {}".format(feature_layer))

        def feature_hook(module, input, output):
            self.feature_blobs.append(output.data.cpu().numpy())

        self.net._modules.get(feature_layer).register_forward_hook(feature_hook)

    def _calc_cam(self, size_upsample):
        # Calculate class activation map.
        feature_conv = self.feature_blobs[0]

        params = list(self.net.parameters())
        weight_softmax = np.squeeze(params[-2].data.numpy())

        logit = self.net(self.matrix)
        prob_and_idx = func.softmax(logit, dim=1).data.squeeze()
        probs, idx = prob_and_idx.sort(0, True)
        probs, idx = probs.numpy, idx.numpy()

        bz, nc, h, w = feature_conv.shape
        cam = weight_softmax[idx[0]].dot(feature_conv.reshape((nc, h*w)))
        cam = cam.reshape(h, w)
        cam = cam - np.min(cam)
        cam_img = cam / np.max(cam)
        cam_img = np.unit8(255 * cam_img)

        self.cam_img = cv2.resize(cam_img, size_upsample)

    def _draw_pic(self, cam_img, alpha=(0.3, 0.5)):
        if isinstance(alpha, float):
            alpha_h = alpha_m = alpha
        elif isinstance(alpha, list) and len(alpha) == 2:
            alpha_h, alpha_m = alpha
        else:
            alpha_h, alpha_m = 0.3, 0.5

        cam_img = cv2.resize(cam_img, (self.matrix.shape))
        self.heatmap = cv2.applyColorMap(cam_img, cv2.COLORMAP_JET)

        return self.heatmap * alpha_h + self.matrix * alpha_m

    def init(self, model_state):
        self.net.load_state_dict(torch.load(model_state))
        self.net.eval()

        self._reg_hook()

        return self

    def show_cam(self, save_path):
        cam_img = self._calc_cam()
        cam_img = self._draw_pic(cam_img)

        cv2.imwrite(save_path, cam_img)
        return self
