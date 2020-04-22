#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author      : Zhenhua Zhang
# Email       : zhenhua.zhang217@gmail.com
# License     : MIT
# Created date: Thu 12 Mar 2020 05:44:52 PM CET
# Last update : Mon 30 Mar 2020 04:32:17 PM CEST

"""A module to train a convolutionary neural network.

```
Example:
    TODO

```
"""

import sys
import math
import glob
import logging
import warnings

warnings.filterwarnings("ignore")

import numpy as np

from sklearn.model_selection import StratifiedKFold
from sklearn.metrics.classification import accuracy_score, precision_score, recall_score

from statsmodels.sandbox.stats.multicomp import multipletests

logger = logging.getLogger("pytorch")
logger.setLevel(logging.ERROR)

logger = logging.getLogger()
fmt = logging.Formatter("| {levelname: ^8} | {asctime} | {name}: {message}", datefmt="%Y-%m-%d %H:%M:%S %p", style="{")
cs_handle = logging.StreamHandler()
cs_handle.setLevel(logging.DEBUG)
cs_handle.setFormatter(fmt)
logger.addHandler(cs_handle)
logger.setLevel(logging.DEBUG)

try:
    import torch
    import torch.nn as nn
    import torch.nn.functional as func
    import torch.optim as optim

    from torch.utils.data import random_split
    from torch.utils.data import Dataset, DataLoader, Subset
except ImportError as ime:
    logger.error("{}".format(ime))


class ASEDataset(Dataset):
    def __init__(self, gene_id, file_pat, element_trans=None, dataset_trans=None):
        """
        Args:
            gene_id   (string): Gene ID (Ensembl gene ID) to train on.
            file_pat  (string): Pattern to find the numpy file.
            element_trans (callable, optional): Optional transfrom to be applied
            on a sample.
            dataset_trans (callable, optional): Optional transfrom to be applied
            on the whole dataset.

        NOTE:
            1. In previous implementation, it's not handy to do multiple-test
            adjustment as the sampel was loaded at getitem. Thus, in current
            implementation, the data will be loaded into memory and then further
            operation will be done.
        """
        self.gene_id = gene_id
        self.file_pat = file_pat
        self.element_trans = element_trans
        self.dataset_trans = dataset_trans

        self.file_path_pool = self._pr_load_file_path()
        self.dataset_pool = self._pr_load_data()

    def __len__(self):
        return len(self.file_path_pool)

    def __getitem__(self, idx):
        sample = self.dataset_pool[idx]

        if self.element_trans:
            sample = self.element_trans(sample)

        return sample

    def _pr_load_file_path(self):
        return glob.glob(self.file_pat, recursive=True)

    def _pr_load_data(self):
        """Load dataset."""
        temp_list = []
        for idx in range(len(self)):
            gene_id = self.gene_id
            if self.file_path_pool:
                file_path = self.file_path_pool[idx]
            else:
                raise ValueError()
            dataset = np.load(file_path, allow_pickle=True).get(self.gene_id)

            if dataset is None:
                logger.error("[E]: Failed to find {} in file: {}".format(gene_id, file_path))
                return None

            matrix = dataset[0].astype(np.float32)
            length = matrix.shape[1]
            max_sampl_num = int(math.sqrt(length / 2)) ** 2 * 2
            if max_sampl_num != length:
                logger.warn("Reshape the input data for the product of length * width has no integer square root solution")
                matrix = matrix[:, :max_sampl_num, :]

            label = dataset[1]
            sample = (matrix, label)
            temp_list.append(sample)

        if self.dataset_trans:
            temp_list = self.dataset_trans(temp_list)

        return tuple(temp_list)

    def _pr_items(self, idx, labels=True):
        pos = 1 if labels else 0
        if idx is None:
            for idx in range(len(self)):
                yield self[idx][pos]
        else:
            yield self[idx][pos]

    def get_labels(self, idx=None):
        return self._pr_items(idx)

    def get_matrix(self, idx=None):
        return self._pr_items(idx, False)


class MultipleTestAdjustMent(object):
    """Adjust p-values for multiple tests.

    This class is a dataset-wide transformer.
    """

    def __init__(self, alpha=0.05, method="fdr_bh"):
        self.alpha = 0.05
        self.method = method

    def __call__(self, dataset: list):
        p_val_list = []
        lable_list = []
        matrix_list = []
        for matrix, (p_val, label) in dataset:
            p_val_list.append(p_val)
            lable_list.append(label)
            matrix_list.append(matrix)
        p_val_list = multipletests(p_val_list, alpha=self.alpha,
                                   method=self.method)[1]

        return tuple(zip(matrix_list, tuple(zip(p_val_list, lable_list))))


class ReshapeMatrixAndPickupLabel(object):
    """Reshape the sequence matrix into a given size.

    This class is an element-wise transformer.
    Args:
        -1, 0, 1 reprsents ASE effect prone to allele A, no ASE effects, ASE
        Effects prone to allele B in the raw ASE quantification results,
        however, both PyTorch and TensorFlow require labels no less than 0.
    """

    def __init__(self, pthd=0.05):
        self.pthd = pthd

    def __call__(self, sample):
        matrix, labels = sample

        n_channel, length, n_type_nucl = matrix.shape
        length = int(math.sqrt(length * n_type_nucl))

        matrix = matrix.reshape((n_channel, length, length))
        if labels[0] < self.pthd:
            label = labels[1] + 1
        else:
            label = 1

        return (matrix, label)


class CNNModel(nn.Module):
    """A built-in CNN model for the package.
    """

    def __init__(self):
        super(CNNModel, self).__init__()

        self.conv01 = nn.Conv2d(1, 8, 3, padding=1)
        self.conv02 = nn.Conv2d(8, 16, 3, padding=1)
        self.conv03 = nn.Conv2d(16, 32, 3, padding=1)
        self.conv04 = nn.Conv2d(32, 64, 3, padding=1)
        self.pool0 = nn.MaxPool2d(3, 1, 1)  # H=128, W=128

        self.conv11 = nn.Conv2d(64, 32, 3, padding=1)
        self.conv12 = nn.Conv2d(32, 16, 3, padding=1)
        self.conv13 = nn.Conv2d(16, 8, 3, padding=1)
        self.conv14 = nn.Conv2d(8, 4, 3, padding=1)
        self.pool1 = nn.MaxPool2d(5, 1)  # H=124, W=124

        self.conv21 = nn.Conv2d(4, 16, 3, padding=1)
        self.conv22 = nn.Conv2d(16, 32, 3, padding=1)
        self.conv23 = nn.Conv2d(32, 128, 3, padding=1)
        self.conv24 = nn.Conv2d(128, 256, 3, padding=1)
        self.pool2 = nn.MaxPool2d(5, 1)  # H=120, W=120

        self.conv31 = nn.Conv2d(256, 128, 3, padding=1)
        self.conv32 = nn.Conv2d(128, 64, 3, padding=1)
        self.conv33 = nn.Conv2d(64, 32, 3, padding=1)
        self.conv34 = nn.Conv2d(32, 16, 3, padding=1)
        self.pool3 = nn.MaxPool2d(5, 1)  # H=116, W=116

        self.conv41 = nn.Conv2d(16, 32, 3, padding=1)
        self.conv42 = nn.Conv2d(32, 64, 3, padding=1)
        self.conv43 = nn.Conv2d(64, 128, 3, padding=1)
        self.conv44 = nn.Conv2d(128, 256, 3, padding=1)
        self.pool4 = nn.MaxPool2d(5, 1)  # H=112, W=112

        self.conv51 = nn.Conv2d(256, 128, 3, padding=1)
        self.conv52 = nn.Conv2d(128, 64, 3, padding=1)
        self.conv53 = nn.Conv2d(64, 32, 3, padding=1)
        self.conv54 = nn.Conv2d(32, 16, 3, padding=1)
        self.pool5 = nn.MaxPool2d(7, 1)  # H=106, W=106

        self.conv61 = nn.Conv2d(16, 32, 3, padding=1)
        self.conv62 = nn.Conv2d(32, 64, 3, padding=1)
        self.conv63 = nn.Conv2d(64, 128, 3, padding=1)
        self.conv64 = nn.Conv2d(128, 256, 3, padding=1)
        self.pool6 = nn.MaxPool2d(9, 1)  # H=98, W=98

        self.conv71 = nn.Conv2d(256, 128, 3, padding=1)
        self.conv72 = nn.Conv2d(128, 64, 3, padding=1)
        self.conv73 = nn.Conv2d(64, 32, 3, padding=1)
        self.conv74 = nn.Conv2d(32, 16, 3, padding=1)
        self.pool7 = nn.MaxPool2d(11, 1)  # H=88, W=88

        self.conv81 = nn.Conv2d(16, 32, 3, padding=1)
        self.conv82 = nn.Conv2d(32, 64, 3, padding=1)
        self.conv83 = nn.Conv2d(64, 128, 3, padding=1)
        self.conv84 = nn.Conv2d(128, 256, 3, padding=1)
        self.pool8 = nn.MaxPool2d(13, 1)  # H=76, W=76

        self.conv91 = nn.Conv2d(256, 128, 3, padding=1)
        self.conv92 = nn.Conv2d(128, 64, 3, padding=1)
        self.conv93 = nn.Conv2d(64, 32, 3, padding=1)
        self.conv94 = nn.Conv2d(32, 4, 3, padding=1)
        self.pool9 = nn.MaxPool2d(13, 1)  # H=64, W=64

        self.flat = nn.Flatten()
        self.fc1 = nn.Linear(64*64*4, 16)
        self.fc2 = nn.Linear(16, 8)
        self.fc3 = nn.Linear(8, 3)

    def forward(self, x):
        x = func.relu(self.conv01(x))
        x = func.relu(self.conv02(x))
        x = func.relu(self.conv03(x))
        x = func.relu(self.conv04(x))
        x = self.pool0(x)

        x = func.relu(self.conv11(x))
        x = func.relu(self.conv12(x))
        x = func.relu(self.conv13(x))
        x = func.relu(self.conv14(x))
        x = self.pool1(x)

        x = func.relu(self.conv21(x))
        x = func.relu(self.conv22(x))
        x = func.relu(self.conv23(x))
        x = func.relu(self.conv24(x))
        x = self.pool2(x)

        x = func.relu(self.conv31(x))
        x = func.relu(self.conv32(x))
        x = func.relu(self.conv33(x))
        x = func.relu(self.conv34(x))
        x = self.pool3(x)

        x = func.relu(self.conv41(x))
        x = func.relu(self.conv42(x))
        x = func.relu(self.conv43(x))
        x = func.relu(self.conv44(x))
        x = self.pool4(x)

        x = func.relu(self.conv51(x))
        x = func.relu(self.conv52(x))
        x = func.relu(self.conv53(x))
        x = func.relu(self.conv54(x))
        x = self.pool5(x)

        x = func.relu(self.conv61(x))
        x = func.relu(self.conv62(x))
        x = func.relu(self.conv63(x))
        x = func.relu(self.conv64(x))
        x = self.pool6(x)

        x = func.relu(self.conv71(x))
        x = func.relu(self.conv72(x))
        x = func.relu(self.conv73(x))
        x = func.relu(self.conv74(x))
        x = self.pool7(x)

        x = func.relu(self.conv81(x))
        x = func.relu(self.conv82(x))
        x = func.relu(self.conv83(x))
        x = func.relu(self.conv84(x))
        x = self.pool8(x)

        x = func.relu(self.conv91(x))
        x = func.relu(self.conv92(x))
        x = func.relu(self.conv93(x))
        x = func.relu(self.conv94(x))
        x = self.pool9(x)

        x = self.flat(x)
        x = func.relu(self.fc1(x))
        x = func.relu(self.fc2(x))

        return func.sigmoid(self.fc3(x))


class DLPFactory:
    """A class to train the neuroal network.
    """

    def __init__(self, net, gene_id=None, file_pat=None, togpu=False):
        self.net = net
        self.togpu = togpu
        self.gene_id = gene_id
        self.file_pat = file_pat
        self.dataset = None
        self.splits = None

    @staticmethod
    def _pr_check_device():
        return "cuda:0" if torch.cuda.is_available() else "cpu"

    def _pr_test(self, model_state: str = None, testloader: DataLoader = None):
        if self.net is None:
            self.net = CNNModel()
            if model_state is not None:
                self.net.load_state_dict(torch.load(model_state))
            else:
                logger.error("The model has NOT trained, require model states")
                sys.exit()

        if testloader is None:
            logger.error("Missing testloader!! Exit...")
            return None

        true_list, pred_list = [], []
        device = self._pr_check_device()
        with torch.no_grad():
            for data in testloader:
                matrix, labels = data[0].to(device), data[1].to(device)
                outputs = self.net(matrix)
                _, predicted = torch.max(outputs.data, 1)

                true_list.extend(labels.to("cpu"))
                pred_list.extend(predicted.to("cpu"))

        precision = precision_score(true_list, pred_list, average="macro") * 100
        recall = recall_score(true_list, pred_list, average="macro") * 100
        accuracy = accuracy_score(true_list, pred_list) * 100

        logger.info("Precision: {:.3}%. Recall: {:.3}%. Accuracy: {:.3}%".format(precision, recall, accuracy))

        return None

    def init(self):
        # Check GPU or CPU
        device = self._pr_check_device()
        self.net.to(device)
        return self

    def load_dataset(self, **kwargs):
        """Load the dataset for train or test.
        """
        gene_id = kwargs.get("gene_id", self.gene_id)
        file_pat = kwargs.get("file_pat", self.file_pat)
        element_trans = kwargs.get("element_trans", ReshapeMatrixAndPickupLabel())
        dataset_trans = kwargs.get("dataset_trans", MultipleTestAdjustMent())

        self.dataset = ASEDataset(file_pat=file_pat, gene_id=gene_id, element_trans=element_trans, dataset_trans=dataset_trans)

        return self

    def k_cv_split(self, n_splits=6, **kwargs):
        """Split the dataset into given number of splits for cross-validation.
        """
        lables = list(self.dataset.get_labels())
        matrix = list(self.dataset.get_matrix())

        skf = StratifiedKFold(n_splits=n_splits)
        self.splits = skf.split(X=matrix, y=lables)

        return self

    def train(self, eps=16, criterion=None, optimizer=None, **kwargs):
        """Train the model."""
        if criterion is None:
            criterion = nn.CrossEntropyLoss()

        if optimizer is None:
            optimizer = optim.Adam(self.net.parameters(), lr=0.0001)

        batch_size = kwargs.get('batch_size', 8)
        shuffle = kwargs.get("shuffle", False)
        num_workers = kwargs.get("num_workers", 4)
        device = self._pr_check_device()
        for cv_idx, (cv_train_idx, cv_test_idx) in enumerate(self.splits):
            trainloader = DataLoader(Subset(self.dataset, cv_train_idx), batch_size=batch_size, shuffle=shuffle, num_workers=num_workers)
            testloader = DataLoader(Subset(self.dataset, cv_test_idx), batch_size=batch_size, shuffle=shuffle, num_workers=num_workers)

            logger.info("Cross validation {}".format(cv_idx))
            for epoch in range(eps):
                running_loss = 0.0
                for _, data in enumerate(trainloader, 0):
                    optimizer.zero_grad()

                    inputs, labels = data[0].to(device), data[1].to(device)
                    outputs = self.net(inputs)
                    loss = criterion(outputs, labels)
                    loss.backward()
                    optimizer.step()

                    running_loss += loss.item()

                if epoch % 10 == 9:
                    logger.info("Epoch: {}, loss: {}".format(epoch+1, running_loss / 10))

                running_loss = 0.0

            self._pr_test(self.net, testloader)

        return self

    def save_model(self, path):
        """Save the trained model, typlically the state dictionary.
        """
        if self.net is None:
            logger.error("The mode is empty, have you ever run the train() method?")
            sys.exit()

        torch.save(self.net.state_dict(), path)

        return self

    def test(self, model_state: str = None, testloader: DataLoader = None):
        """Test for new dataset.
        """
        self._pr_test(model_state, testloader)

        return self

    def predict(self, matrix=None):
        """Apply the model on given dataset.
        """
        logger.warnings("Note implemented yet!")

        return self
