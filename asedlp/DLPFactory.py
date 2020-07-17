#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author      : Zhenhua Zhang
# Email       : zhenhua.zhang217@gmail.com
# License     : MIT
# Created date: Thu 12 Mar 2020 05:44:52 PM CET

"""A module to train a convolutionary neural network."""

import sys

import numpy as nmp
import seaborn as sbn
import matplotlib.pyplot as plt
from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.metrics import accuracy_score, precision_score, recall_score, roc_auc_score, roc_curve

import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as func
from torch.utils.data import DataLoader, Subset
from torch.utils.tensorboard import SummaryWriter

from ASEDataset import ASEDataset, ReshapeMatrixAndPickupLabel, MultipleTestAdjustMent
from zutils import logger


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
    def _check_device():
        return "cuda:0" if torch.cuda.is_available() else "cpu"

    @staticmethod
    def _eval_metric(true_label, pred_label, pred_score, class_label=[0, 1, 2], percent=True):
        scale = 100 if percent else 1

        pred_score = func.softmax(pred_score, dim=1).data.numpy()
        fpr_list, tpr_list = [], []
        for idx, label in enumerate(class_label):
            true_label_binary = [1 if x == label else -1 for x in true_label]
            pred_score_ = pred_score[:, idx]
            fpr, tpr, _ = roc_curve(true_label_binary, pred_score_)

            fpr_list.append(fpr)
            tpr_list.append(tpr)

        rcac = roc_auc_score(true_label, pred_score, average="macro", multi_class="ovr") * scale
        prec = precision_score(true_label, pred_label, average="macro") * scale
        rcll = recall_score(true_label, pred_label, average="macro") * scale
        accu = accuracy_score(true_label, pred_label) * scale

        return fpr_list, tpr_list, rcac, prec, accu, rcll

    def _test(self, model_state: str = None, testloader: DataLoader = None,
                 criterion=None):
        if self.net is None:
            if model_state is not None:
                self.net.load_state_dict(torch.load(model_state))
            else:
                logger.error("The model has NOT trained, require model states")
                sys.exit()

        if testloader is None:
            logger.error("Missing testloader!! Exit...")
            return None, None, None, None, None, None

        if criterion is None:
            criterion = nn.CrossEntropyLoss()

        true_label_list, pred_label_list, pred_score_list = [], [], None
        device = self._check_device()
        self.net.to(device)
        test_loss = 0.0
        with torch.no_grad():
            for data in testloader:
                inputs, true_label = data[0].to(device), data[1].to(device)
                outputs = self.net(inputs)

                _, pred_label = torch.max(outputs.data, 1)
                pred_label_list.extend(pred_label.to("cpu"))
                true_label_list.extend(true_label.to("cpu"))
                
                if pred_score_list is None:
                    pred_score_list = outputs.data.to("cpu")
                else:
                    pred_score_list = torch.cat((pred_score_list, outputs.data.to("cpu")))

                test_loss += criterion(outputs, true_label).item()

        return true_label_list, pred_label_list, pred_score_list, test_loss

    def _add_roc_curve(self, fpr_list, tpr_list, tag):
        for class_idx in range(len(fpr_list)):
            fpr = fpr_list[class_idx] * 100
            tpr = tpr_list[class_idx] * 100
            for idx, x in enumerate(fpr):
                self.writer.add_scalars("ROC_AUC_curve/{}".format(tag), {str(class_idx): tpr[idx]}, x)

    def _train(self, eps, optimizer, device, criterion, splits, batch_size, shuffle, num_workers, log_per_n_epoch=5):
        """The exact implementation of train.
        """
        logger.info("| CV | Operation | Epoch | Accuracy | Precision | Recall | ROC_AUC | Loss")
        for cv_idx, (train_idx, test_idx) in enumerate(splits):
            trainloader = DataLoader(Subset(self.dataset, train_idx), batch_size=batch_size, shuffle=shuffle, num_workers=num_workers)
            testloader = DataLoader(Subset(self.dataset, test_idx), batch_size=batch_size, shuffle=shuffle, num_workers=num_workers)

            cv_loss_list = []
            for epoch in range(eps):
                running_loss = 0.0
                true_label_list, pred_label_list, pred_score_list = [], [], None
                for _, data in enumerate(trainloader):
                    optimizer.zero_grad()

                    inputs, true_label = data[0].to(device), data[1].to(device)
                    outputs = self.net(inputs)
                    train_loss = criterion(outputs, true_label)
                    train_loss.backward()
                    optimizer.step()

                    _, pred_label = torch.max(outputs.data, 1)
                    pred_label_list.extend(pred_label.to("cpu"))
                    true_label_list.extend(true_label.to("cpu"))

                    if pred_score_list is None:
                        pred_score_list = outputs.data.to("cpu")
                    else:
                        pred_score_list = torch.cat((pred_score_list, outputs.data.to("cpu")))

                    running_loss += train_loss.item()

                if epoch == 0 and cv_idx == 0: # 确保只添加一次graph
                    self.writer.add_graph(self.net, inputs)

                if (epoch + 1) % log_per_n_epoch == 0 or epoch < 10:
                    tn_fpr, tn_tpr, tn_roc_auc, tn_prec, tn_accu, tn_rcll = self._eval_metric(true_label_list, pred_label_list, pred_score_list)
                    self.writer.add_scalar("ROC_AUC_score/Train", tn_roc_auc, epoch)
                    self.writer.add_scalar("Precision/Train", tn_prec, epoch)
                    self.writer.add_scalar("Accuracy/Train", tn_accu, epoch)
                    self.writer.add_scalar("Recall/Train", tn_rcll, epoch)
                    self.writer.add_scalar("Loss/Train", running_loss, epoch)
                    logger.info("| {} | Train | {: ^3} | {: ^5.2f}% | {: ^5.2f}% | {: ^5.2f}% | {: ^5.2f}% | {:.4f}".format(cv_idx + 1, epoch+1, tn_accu, tn_prec, tn_rcll, tn_roc_auc, running_loss))

                    true_label_list, pred_label_list, pred_score_list, tt_loss = self._test(self.net, testloader)
                    tt_fpr, tt_tpr, tt_roc_auc, tt_prec, tt_accu, tt_rcll = self._eval_metric(true_label_list, pred_label_list, pred_score_list)
                    self.writer.add_scalar("ROC_AUC_score/Test", tt_roc_auc, epoch)
                    self.writer.add_scalar("Precision/Test", tt_prec, epoch)
                    self.writer.add_scalar("Accuracy/Test", tt_accu, epoch)
                    self.writer.add_scalar("Recall/Test", tt_rcll, epoch)
                    self.writer.add_scalar("Loss/Test", tt_loss, epoch)
                    logger.info("| {} | Test  | {: ^3} | {: ^5.2f}% | {: ^5.2f}% | {: ^5.2f}% | {: ^5.2f}% | {:.4f}".format(cv_idx+1, epoch+1, tt_accu, tt_prec, tt_rcll, tt_roc_auc, tt_loss))

                cv_loss_list.append(running_loss)

        self._add_roc_curve(tn_fpr, tn_tpr, "Train")
        self._add_roc_curve(tt_fpr, tt_tpr, "Test")
        self.train_loss_list.append(cv_loss_list)

    def init(self):
        """Init."""
        # Check GPU or CPU
        self.writer = SummaryWriter(self.logging_path)
        self.net.to(self._check_device())
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
        labels = list(self.dataset.get_labels())
        matrix = list(self.dataset.get_matrix())

        skf = StratifiedKFold(n_splits=n_splits)
        self.cv_splits = skf.split(X=matrix, y=labels)

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
        device = self._check_device()

        if self.cv_splits is None and self.tt_splits is None:
            logger.error("No train and test splits were found.")
            return self

        if self.tt_splits is None:
            splits = self.cv_splits
        else:
            splits = [self.tt_splits]

        log_per_n_epoch = self.log_per_n_epoch
        self._train(eps, optimizer, device, criterion, splits, batch_size,
                       shuffle, num_workers, log_per_n_epoch=log_per_n_epoch)

        self.writer.close()
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
        self._test(model_state, testloader)

        return self

    def predict(self, matrix=None):
        """Apply the model on given dataset.
        """
        logger.warnings("Note implemented yet!")

        return self
