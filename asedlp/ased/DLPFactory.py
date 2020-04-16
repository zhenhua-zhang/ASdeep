#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author      : Zhenhua Zhang
# Email       : zhenhua.zhang217@gmail.com
# License     : MIT
# Created date: Thu 12 Mar 2020 05:44:52 PM CET
# Last update : Mon 30 Mar 2020 04:32:17 PM CEST

"""A module to train a convolutionary neural network."""
import sys
import math
import glob
import logging
import warnings

import numpy as np
import sklearn.model_selection as skms

from statsmodels.sandbox.stats.multicomp import multipletests

warnings.filterwarnings("ignore")

logger = logging.getLogger("pytorch")
logger.setLevel(logging.ERROR)

logger = logging.getLogger()
fmt = logging.Formatter("| {levelname: ^8} | {asctime} | {name}: {message}", datefmt="%Y-%m-%d %H:%M:%S %p", style="{")
cs_handle = logging.StreamHandler()
cs_handle.setLevel(logging.INFO)
cs_handle.setFormatter(fmt)
logger.addHandler(cs_handle)
logger.setLevel(logging.INFO)

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
            element_trans (callable, optional): Optional transfrom to be applied on
            a sample.
            dataset_trans (callable, optional): Optional transfrom to be applied on
            the whole dataset.

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
    """Adjust p-values for multiple tests."""

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
    def __init__(self):
        super(CNNModel, self).__init__()
        self.pool = nn.MaxPool2d(2, 2)
        self.conv1 = nn.Conv2d(1, 4, 3)
        self.conv2 = nn.Conv2d(4, 8, 3)
        self.conv3 = nn.Conv2d(8, 16, 3)
        self.conv4 = nn.Conv2d(16, 32, 3)
        self.conv5 = nn.Conv2d(32, 16, 3)
        self.flat = nn.Flatten()
        self.fc1 = nn.Linear(2*2*16, 1024)
        self.fc2 = nn.Linear(1024, 512)
        self.fc3 = nn.Linear(512, 256)
        self.fc4 = nn.Linear(256, 128)
        self.fc5 = nn.Linear(128, 64)
        self.fc6 = nn.Linear(64, 3)

    def forward(self, x):
        x = self.pool(func.relu(self.conv1(x)))
        x = self.pool(func.relu(self.conv2(x)))
        x = self.pool(func.relu(self.conv3(x)))
        x = self.pool(func.relu(self.conv4(x)))
        x = self.pool(func.relu(self.conv5(x)))
        x = self.flat(x)
        x = func.relu(self.fc1(x))
        x = func.relu(self.fc2(x))
        x = func.relu(self.fc3(x))
        x = func.relu(self.fc4(x))
        x = func.relu(self.fc5(x))
        return func.sigmoid(self.fc6(x))


class DLPFactory:
    def __init__(self, net, gene_id=None, file_pat=None, signature="train", togpu=False):
        self.net = net
        self.togpu = togpu
        self.gene_id = gene_id
        self.file_pat = file_pat
        self.signature = signature
        self.dataset = None
        self.splits = None

    def _pr_check_device(self):
        return "cuda:0" if torch.cuda.is_available() else "cpu"

    def _pr_test(self, model_state: str = None, testloader: DataLoader = None):
        if self.net is None:
            self.net = CNNModel()
            if model_state is not None:
                self.net.load_state_dict(torch.load(model_state))
            else:
                logger.error("The model has NOT trained, please supply the model state")
                sys.exit()

        correct = 0
        total = 0

        if testloader is None:
            logger.error("Missing testloader!! Exit...")
            return self

        with torch.no_grad():
            for data in testloader:
                matrix, labels = data
                outputs = self.net(matrix)

                _, predicted = torch.max(outputs.data, 1)
                total += labels.size(0)
                correct += (predicted == labels).sum().item()

        logger.info("Accuracy of the network on the test matrix: {} %".format(100 * correct / total))

    def init(self):
        # Check GPU or CPU
        device = self._pr_check_device()
        self.net.to(device)

        return self

    def load_dataset(self, train_prop=0.8, **kwargs):
        """Load the dataset for train or test."""
        # Create a dataset
        gene_id = kwargs.get("gene_id", self.gene_id)
        file_pat = kwargs.get("file_pat", self.file_pat)
        element_trans = kwargs.get("element_trans", ReshapeMatrixAndPickupLabel())
        dataset_trans = kwargs.get("dataset_trans", MultipleTestAdjustMent())

        self.dataset = ASEDataset(file_pat=file_pat, gene_id=gene_id,
                                  element_trans=element_trans,
                                  dataset_trans=dataset_trans)

        return self

    def k_cv_split(self, n_splits=6, **kwargs):
        """Split the dataset into given number of splits for cross-validation.
        """
        lables = list(self.dataset.get_labels())
        index = range(len(self.dataset))

        skf = skms.StratifiedKFold(n_splits=n_splits)
        splits = skf.get_n_splits(X=index, y=lables)
        self.splits = splits

        return self

    def train(self, eps=16, criterion=None, optimizer=None, **kwargs):
        """Train the model."""
        if criterion is None:
            criterion = nn.CrossEntropyLoss()
        if optimizer is None:
            optimizer = optim.Adam(self.net.parameters(), lr=0.001)

        device = self._pr_check_device()

        batch_size = kwargs.get('batch_size', 8)
        shuffle = kwargs.get("shuffle", False)
        num_workers = kwargs.get("num_workers", 4)

        for cv_idx, (cv_train_idx, cv_test_idx) in enumerate(self.splits):
            trainloader = DataLoader(
                Subset(self.dataset, cv_train_idx), batch_size=batch_size,
                shuffle=shuffle, num_workers=num_workers)  # Train dataset
            testloader = DataLoader(
                Subset(self.dataset, cv_test_idx), batch_size=batch_size,
                shuffle=shuffle, num_workers=num_workers)   # Test dataset

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

            self._pr_test(testloader)

        return self

    def save_model(self, path):
        if self.net is None:
            logger.error("The mode is empty, have you ever run the train() method?")
            sys.exit()

        torch.save(self.net.state_dict(), path)
        return self

    def test(self, model_state: str = None, testloader: DataLoader = None):
        self._pr_test(model_state, testloader)

        return self
