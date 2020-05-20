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

import logging
import math
import sys
import warnings

import matplotlib.pyplot as plt
import numpy as np
from statsmodels.sandbox.stats.multicomp import multipletests

import seaborn as sbn
from sklearn.metrics.classification import accuracy_score, precision_score, recall_score
from sklearn.model_selection import StratifiedKFold, train_test_split

warnings.filterwarnings("ignore")

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
    from torch.utils.data import Dataset, DataLoader, Subset
    from torch.utils.tensorboard import SummaryWriter
except ImportError as ime:
    logger.error("{}".format(ime))


class ASEDataset(Dataset):
    def __init__(self, gene_id, file_path_pool, element_trans=None, dataset_trans=None):
        """
        Args:
            gene_id   (string): Gene ID (Ensembl gene ID) to train on.
            file_path_pool  (string): Pattern to find the numpy file.
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
        self.element_trans = element_trans
        self.dataset_trans = dataset_trans
        self.file_path_pool = file_path_pool


        self.dataset_pool = self._pr_load_data()

    def __len__(self):
        return len(self.file_path_pool)

    def __getitem__(self, idx):
        sample = self.dataset_pool[idx]

        if self.element_trans:
            sample = self.element_trans(sample)

        return sample

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

    The current implementation is VGG16.
    """

    def __init__(self, input_size=128, dropout_p=0.5):
        super(CNNModel, self).__init__()

        self.pool = nn.MaxPool2d(2, 2)
        self.conv01 = nn.Conv2d(1, 64, 3, stride=1, padding=1)
        self.conv02 = nn.Conv2d(64, 64, 3, stride=1, padding=1)

        self.conv03 = nn.Conv2d(64, 128, 3, stride=1, padding=1)
        self.conv04 = nn.Conv2d(128, 128, 3, stride=1, padding=1)

        self.conv05 = nn.Conv2d(128, 256, 3, stride=1, padding=1)
        self.conv06 = nn.Conv2d(256, 256, 3, stride=1, padding=1)
        self.conv07 = nn.Conv2d(256, 256, 3, stride=1, padding=1)

        self.conv08 = nn.Conv2d(256, 512, 3, stride=1, padding=1)
        self.conv09 = nn.Conv2d(512, 512, 3, stride=1, padding=1)
        self.conv10 = nn.Conv2d(512, 512, 3, stride=1, padding=1)

        self.conv11 = nn.Conv2d(512, 512, 3, stride=1, padding=1)
        self.conv12 = nn.Conv2d(512, 512, 3, stride=1, padding=1)
        self.conv13 = nn.Conv2d(512, 512, 3, stride=1, padding=1)

        self.flat = nn.Flatten()
        self.drop = nn.Dropout(p=dropout_p)

        self.fc1 = nn.Linear(int(512 * (input_size / 32) ** 2), 4096)
        self.fc2 = nn.Linear(4096, 4096)
        self.fc3 = nn.Linear(4096, 3)

    def forward(self, x):
        """forward propogation."""
        x = func.relu(self.conv01(x))
        x = func.relu(self.conv02(x))
        x = self.pool(x)

        x = func.relu(self.conv02(x))
        x = func.relu(self.conv03(x))
        x = self.pool(x)

        x = func.relu(self.conv05(x))
        x = func.relu(self.conv06(x))
        x = func.relu(self.conv07(x))
        x = self.pool(x)

        x = func.relu(self.conv08(x))
        x = func.relu(self.conv09(x))
        x = func.relu(self.conv10(x))
        x = self.pool(x)

        x = func.relu(self.conv11(x))
        x = func.relu(self.conv12(x))
        x = func.relu(self.conv13(x))
        x = self.pool(x)

        x = self.flat(x)
        x = func.relu(self.fc1(x))
        x = self.drop(x)
        x = func.relu(self.fc2(x))
        x = self.drop(x)

        return func.sigmoid(self.fc3(x))


class DLPFactory:
    """A class to train the neuroal network.
    """

    def __init__(self, net, gene_id=None, file_path_pool=None, togpu=False,
                 logging_path=None, log_per_n_epoch=5):
        self.net = net                         # NN model will be used
        self.togpu = togpu                     # whether ship data and model to GPU
        self.gene_id = gene_id                 # The gene ID on which will train the model
        self.logging_path = logging_path       # Path for the TensorBoard logs
        self.file_path_pool = file_path_pool   # Files on which train and test
        self.log_per_n_epoch = log_per_n_epoch # Logging point per N epoches

        self.dataset = None    # DataSet() for the training
        self.cv_splits = None  # cross validation splits
        self.tt_splits = None  # train:test splits

        self.train_loss_list = []  # Loss

        self.writer = None # SummaryWriter for TensorBoard

    @staticmethod
    def _pr_check_device():
        return "cuda:0" if torch.cuda.is_available() else "cpu"

    def _pr_train(self, eps, optimizer, device, criterion, splits, batch_size, shuffle, num_workers, log_per_n_epoch=5):
        """The exact implementation of train.
        """
        for cv_idx, (train_idx, test_idx) in enumerate(splits):
            trainloader = DataLoader(Subset(self.dataset, train_idx), batch_size=batch_size, shuffle=shuffle, num_workers=num_workers)
            testloader = DataLoader(Subset(self.dataset, test_idx), batch_size=batch_size, shuffle=shuffle, num_workers=num_workers)

            logger.info("Cross validation {}".format(cv_idx + 1))
            cv_loss_list = []
            for epoch in range(eps):
                running_loss = 0.0
                for _, data in enumerate(trainloader):
                    optimizer.zero_grad()

                    inputs, labels = data[0].to(device), data[1].to(device)
                    outputs = self.net(inputs)
                    train_loss = criterion(outputs, labels)
                    train_loss.backward()
                    optimizer.step()

                    if epoch == 0 and cv_idx == 0:
                        self.writer.add_graph(self.net, inputs)

                    running_loss += train_loss.item()

                if (epoch + 1) % log_per_n_epoch == 0 or epoch < 10:
                    _, predicted = torch.max(outputs.data, 1)
                    train_accu = accuracy_score(labels.to("cpu"), predicted.to("cpu")) * 100
                    train_prec = precision_score(labels.to("cpu"), predicted.to("cpu"), average="macro") * 100
                    train_rcll = recall_score(labels.to("cpu"), predicted.to("cpu"), average="macro") * 100
                    logger.info("TRAIN: epoch {}, accuracy {:.3}%, precision {:.3}%, recall {:.3}%, loss {:.3}".format(epoch+1, train_accu, train_prec, train_rcll, running_loss))

                    self.writer.add_scalar("Accuracy/Train", train_accu, epoch)
                    self.writer.add_scalar("Precision/Train", train_prec, epoch)
                    self.writer.add_scalar("Recall/Train", train_rcll, epoch)
                    self.writer.add_scalar("Loss/Train", running_loss, epoch)

                    test_accu, test_prec, test_rcll, test_loss = self._pr_test(self.net, testloader)
                    logger.info("TEST: epoch {}, accuracy {:.3}%, precision {:.3}%, recall {:.3}%, loss {:.3}".format(epoch+1, test_accu, test_prec, test_rcll, test_loss))

                    self.writer.add_scalar("Accuracy/Test", test_accu, epoch)
                    self.writer.add_scalar("Precision/Test", test_prec, epoch)
                    self.writer.add_scalar("Recall/Test", test_rcll, epoch)
                    self.writer.add_scalar("Loss/Test", test_loss, epoch)

                cv_loss_list.append(running_loss)

        self.train_loss_list.append(cv_loss_list)

    def _pr_test(self, model_state: str = None, testloader: DataLoader = None,
                 criterion=None):
        if self.net is None:
            self.net = CNNModel()
            if model_state is not None:
                self.net.load_state_dict(torch.load(model_state))
            else:
                logger.error("The model has NOT trained, require model states")
                sys.exit()

        if testloader is None:
            logger.error("Missing testloader!! Exit...")
            return None, None, None

        if criterion is None:
            criterion = nn.CrossEntropyLoss()

        true_list, pred_list = [], []
        device = self._pr_check_device()
        self.net.to(device)
        test_loss = 0.0
        with torch.no_grad():
            for data in testloader:
                inputs, labels = data[0].to(device), data[1].to(device)
                outputs = self.net(inputs)
                _, predicted = torch.max(outputs.data, 1)

                true_list.extend(labels.to("cpu"))
                pred_list.extend(predicted.to("cpu"))

                test_loss += criterion(outputs, labels).item()

        precision = precision_score(true_list, pred_list, average="macro") * 100
        recall = recall_score(true_list, pred_list, average="macro") * 100
        accuracy = accuracy_score(true_list, pred_list) * 100

        # logger.info("Precision: {:.3}%. Recall: {:.3}%. Accuracy: {:.3}%. Test loss: {:.3}".format(precision, recall, accuracy, test_loss))

        return accuracy, precision, recall, test_loss

    def init(self):
        """Init."""
        # Check GPU or CPU
        self.writer = SummaryWriter(self.logging_path)
        self.net.to(self._pr_check_device())
        return self

    def load_dataset(self, **kwargs):
        """Load the dataset for train or test.
        """
        gene_id = kwargs.get("gene_id", self.gene_id)
        file_path_pool = kwargs.get("file_path_pool", self.file_path_pool)
        element_trans = kwargs.get("element_trans", ReshapeMatrixAndPickupLabel())
        dataset_trans = kwargs.get("dataset_trans", MultipleTestAdjustMent())

        self.dataset = ASEDataset(file_path_pool=file_path_pool, gene_id=gene_id, element_trans=element_trans, dataset_trans=dataset_trans)

        return self

    def k_cv_split(self, n_splits=10):
        """Split the dataset into given number of splits for cross-validation.

        TODO: 目标是用交叉验证，但实际效果是增加了n_splits倍的epoch。需要改进算法。
        """
        lables = list(self.dataset.get_labels())
        matrix = list(self.dataset.get_matrix())

        skf = StratifiedKFold(n_splits=n_splits)
        self.cv_splits = skf.split(X=matrix, y=lables)

        return self

    def train_test_split(self, train_size=0.7):
        """Split the dataset into train and test dataset."""
        labels = list(self.dataset.get_labels())
        label_idx = range(len(labels))

        self.tt_splits = train_test_split(label_idx, train_size=train_size, stratify=labels)

        return self

    def train(self, eps=20, criterion=None, optimizer=None, learning_rate=1e-5, **kwargs):
        """Train the model.

        Parameters
        ----------
        eps (int, optional, 20): Epoches.
        criterion (object, optional, None): Method to calculate loss.
        optimizer (object, optional, None): Method to optimize the model.
        learning_rate (float, optional, 1e-5): Learning rate. Works only when `optimizer` is None.
        **kwargs (Any, optinal): Other keyword arguments.
        """
        if criterion is None:
            criterion = nn.CrossEntropyLoss()

        if optimizer is None:
            optimizer = optim.Adam(self.net.parameters(), lr=learning_rate)

        batch_size = kwargs.get('batch_size', 32)
        shuffle = kwargs.get("shuffle", False)
        num_workers = kwargs.get("num_workers", 4)
        device = self._pr_check_device()

        if self.cv_splits is None and self.tt_splits is None:
            logger.error("No train and test splits were found.")
            return self

        if self.tt_splits is None:
            splits = self.cv_splits
        else:
            splits = [self.tt_splits]

        log_per_n_epoch = self.log_per_n_epoch
        self._pr_train(eps, optimizer, device, criterion, splits, batch_size,
                       shuffle, num_workers, log_per_n_epoch=log_per_n_epoch)

        return self

    def loss_curve(self, loss_curve_path=None, svfmt="png", **kwargs):
        """Save the loss curve."""
        train_loss_list = self.train_loss_list
        n_cv = len(train_loss_list)

        fig, axes = plt.subplots(n_cv)
        if n_cv == 1:
            loss_axe_pair = [[train_loss_list[0], axes]]
        elif n_cv > 1:
            loss_axe_pair = zip(train_loss_list, axes)
        else:
            logger.error("The loss list is empty!!")
            return self

        for loss, ax in loss_axe_pair:
            epoch_list = list(range(len(loss)))
            sbn.lineplot(x=epoch_list, y=loss, ax=ax)
            ax.set_title("Loss per epoch")
            ax.set_xlabel("Epoch")
            ax.set_ylabel("Loss")

        if loss_curve_path.endswith("/"):
            loss_curve_path = "{path}loss_curve.{fmt}".format(path=loss_curve_path, fmt=svfmt)
        else:
            loss_curve_path = "{path}.loss_curve.{fmt}".format(path=loss_curve_path, fmt=svfmt)

        fig.savefig(loss_curve_path, **kwargs)

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
