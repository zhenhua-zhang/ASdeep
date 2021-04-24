#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''An implementation of Class Activation Mapping (CAM) based on 

Reference: https://github.com/zhoubolei/CAM/blob/master/pytorch_CAM.py
'''

import logging

import torch
import numpy as np
import torch.nn as nn
import matplotlib.pyplot as plt

from PIL import Image
from torchvision import models
from torch.autograd import Variable
from torch.nn import functional as func

from .dataset import ASEDataset, SeqToHelbertAndMakeLabel, MultipleTestAdjustment


class Predictor:
    def __init__(self, net=None, gene_id=None, file_path=None, feature_layer='layer4'):
        self.net = net
        self.gene_id = gene_id
        self.file_path = file_path
        self.feature_layer = feature_layer

        self._feature_blobs = []
        self._matrix = None
        self._hcurve = None
        self._pred_prob = None

        self._true_label = None
        self._pred_label = None

        self._dataset = None
        self._cam_img = None
        self._seq_img = []
        self._weight_softmax = None

    @staticmethod
    def _check_device():
        return 'cuda:0' if torch.cuda.is_available() else 'cpu'

    def _check_matrix_idx(self, matrix_idx):
        if matrix_idx is None:
            return range(len(self._dataset))
        elif isinstance(matrix_idx, int):
            return [matrix_idx]
        elif isinstance(matrix_idx, (list, tuple)) or hasattr(matrix_idx, '__next__'):
            return matrix_idx
        else:
            raise TypeError('Unsupported matrix_idx: {}.'.format(type(matrix_idx)))

    def _reg_hook(self, hook=None, flayer='layer4'):
        def feature_hook(module, input, output):
            self._feature_blobs.append(output.data.cpu().numpy())

        if flayer not in self.net._modules:
            raise KeyError('Unknown feature layer: {}'.format(flayer))

        if hook is None:
            hook = feature_hook

        self.net._modules.get(flayer).register_forward_hook(hook)

    def _load_mtrx(self, matrix_idx):
        self._hcurve, self._true_label = self._dataset[matrix_idx]
        self._matrix = self._hcurve.get_hcurve(onehot=False)

    def _predict(self):
        # Predict the given dataset.
        matrix = Variable(torch.Tensor(np.expand_dims(self._matrix, 0)))
        logit = self.net(matrix)
        probs, idx = func.softmax(logit, dim=1).data.squeeze().sort(0, True)
        self._pred_prob, self._pred_label = probs.numpy()[0], idx.numpy()[0]

    def _calc_cam(self):
        # Calculate class activation map.
        _, nc, h, w = self._feature_blobs[-1].shape
        cam = self._weight_softmax[self._pred_label].dot(self._feature_blobs[-1].reshape((nc, h*w)))
        cam = cam.reshape(h, w)
        cam = np.array(Image.fromarray(cam).resize(tuple(self._matrix.shape[-2:])))
        cam = cam - np.min(cam)
        cam = cam / np.max(cam)

        self._cam_img = np.int16(cam * 255)

    def _cam_to_seq(self):
        for i, j in zip(self._hcurve.x_pool, self._hcurve.y_pool):
            i, j = int(i), int(j)
            self._seq_img.append(self._cam_img[i][j])

        self._seq_img = np.array(self._seq_img)

    def _onehot_to_text(self):
        # Placeholder to implement function convert one-hot encoded sequence
        # into nucleotide sequence.
        self._seq_img

    def _draw_cam_along_seq(self, save_path, title):
        fig, ax = plt.subplots()
        height = self._seq_img.shape[0]
        average_cam = self._seq_img    # 255 is the scale factor

        ax.plot(range(height), average_cam)
        ax.set_title(title)
        ax.set_xlabel('Genomic coordination')
        ax.set_ylabel('CAM score (normalized)')
        fig.savefig(save_path)
        plt.close(fig)

    def init(self, model_state):
        device = self._check_device()

        if self.net is None: # Load model
            self.net = models.resnext50_32x4d()
            self.net.conv1 = nn.Conv2d(1, 64, 7, 2, 3, bias=False) # in_channels 1, out_channels 64, kernel_size 7, stride 2, padding 3
            self.net.fc = nn.Linear(2048, 3)

        self.net.load_state_dict(torch.load(model_state, map_location=device))
        self.net.eval()

        # Fetch parameters
        params = list(self.net.parameters())
        self._weight_softmax = np.squeeze(params[-2].data.numpy())

        self._reg_hook(flayer=self.feature_layer)

        return self

    def load_dataset(self, **kwargs):
        '''Load dataset.'''
        gene_id = kwargs.get('gene_id', self.gene_id)
        file_path_pool = kwargs.get('file_path_pool', self.file_path)
        element_trans = kwargs.get('element_trans', SeqToHelbertAndMakeLabel(onehot=False, matrix=False))
        dataset_trans = kwargs.get('dataset_trans', MultipleTestAdjustment())

        self._dataset = ASEDataset(file_path_pool=file_path_pool, gene_id=gene_id, element_trans=element_trans, dataset_trans=dataset_trans)

        return self

    def predict(self, save_path, matrix_idx=None, show_cam=False):
        matrix_idx = self._check_matrix_idx(matrix_idx)
        lable_dict = {0: 'ASEtoA', 1: 'NonASE', 2: 'ASEtoB'}

        logging.info('TrueLabel | PredLabel | PredProb')
        for idx in matrix_idx:
            self._load_mtrx(idx)
            self._predict()
            logging.info('{} | {} | {}'.format(self._true_label, self._pred_label, self._pred_prob))

            if show_cam:
                self._calc_cam()
                self._cam_to_seq()

                _save_name = save_path + '.{}_{}_'.format(idx, self._true_label)
                self._hcurve.hcurve_to_img(_save_name, overlay=self._cam_img)

                true_label, pred_label = lable_dict[self._true_label], lable_dict[self._pred_label]
                pic_title = 'sample-{:0>3}: {}(true label), {}(pred label)'.format(idx, true_label, pred_label)
                _save_name = save_path + '.{}_{}_{}_cam_along_seq.pdf'.format(idx, self._true_label, self._pred_label)
                self._draw_cam_along_seq(_save_name, pic_title)

        return self


if __name__ == "__main__":
    logging.warning("This module should not be executed directly.")
