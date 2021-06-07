#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author        : Zhenhua Zhang
# Email         : zhenhua.zhang217@gmail.com
# License     : MIT
# Created date: Thu 12 Mar 2020 05:44:52 PM CET

'''A module to train a convolutionary neural network.'''

import sys
import logging

import numpy as np
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve

import torch
import torch.nn as nn
import torch.optim as optim
# import torch.autograd as autograd
import torch.nn.functional as func
from torch.utils.data import DataLoader, Subset
from torch.utils.tensorboard import SummaryWriter

from .dataset import ASEDataset, SeqToHilbertAndMakeLabel, MultipleTestAdjustment


class Trainer:
    '''A class to train the neuroal network.'''
    def __init__(self, net, gene_id=None, file_path_pool=None, logging_path=None,
                 log_per_n_epoch=5, num_workers=12):
        self.net = net                         # NN model will be used
        self.gene_id = gene_id                 # The gene ID on which will train the model
        self.num_workers = num_workers         # Number of workers to load data
        self.logging_path = logging_path       # Path for the TensorBoard logs
        self.file_path_pool = file_path_pool   # Files on which train and test
        self.log_per_n_epoch = log_per_n_epoch # Logging point per N epoches

        self.tt_splits = None # train:test splits

        self.train_loss_list = [] # Loss

        self.writer: SummaryWriter
        self.dataset: ASEDataset

    @staticmethod
    def _check_device():
        return 'cuda:0' if torch.cuda.is_available() else 'cpu'

    @staticmethod
    def _eval_metric(true_label, pred_label, pred_score, class_label=[0, 1, 2], percent=True):
        scale = 100 if percent else 1

        pred_score = func.softmax(pred_score, dim=1).data.numpy()
        fpr_list, tpr_list = [], []
        for idx, label in enumerate(class_label):
            true_label_binary = [1 if x == label else -1 for x in true_label]
            _pred_score = pred_score[:, idx]
            fpr, tpr, _ = roc_curve(true_label_binary, _pred_score)

            fpr_list.append(fpr)
            tpr_list.append(tpr)

        onehot_true_labels = np.zeros((len(true_label), 3))
        for item_idx, class_idx in enumerate(true_label):
            onehot_true_labels[item_idx][class_idx] = 1

        try:
            rcac = roc_auc_score(onehot_true_labels, pred_score, average='macro', multi_class='ovr') * scale
        except ValueError as e:
            rcac = 0

        prec = precision_score(true_label, pred_label, average='macro', zero_division=1) * scale
        rcll = recall_score(true_label, pred_label, average='macro', zero_division=1) * scale
        accu = accuracy_score(true_label, pred_label) * scale

        return fpr_list, tpr_list, rcac, prec, accu, rcll

    def _test(self, model_state: str = None, testloader: DataLoader = None, criterion=None):
        if self.net is None:
            if model_state is not None:
                self.net.load_state_dict(torch.load(model_state))
            else:
                logging.error('The model has NOT trained, require model states')
                sys.exit()

        if testloader is None:
            logging.error('Missing testloader!! Exit...')
            return None, None, None, None

        if criterion is None:
            criterion = nn.CrossEntropyLoss()

        true_label_list, pred_label_list, pred_score_list = [], [], None
        device = self._check_device()
        self.net.to(device)
        test_loss = 0.0
        with torch.no_grad(): # The model should be trained.
            for data in testloader:
                inputs_type = torch.FloatTensor if device == 'cpu' else torch.cuda.FloatTensor

                inputs = data[0].type(inputs_type).to(device)
                true_label = data[1][1].to(device)

                outputs = self.net(inputs)

                _, pred_label = torch.max(outputs.data, 1)
                pred_label_list.extend(pred_label.to('cpu'))
                true_label_list.extend(true_label.to('cpu'))
                
                if pred_score_list is None:
                    pred_score_list = outputs.data.to('cpu')
                else:
                    pred_score_list = torch.cat((pred_score_list, outputs.data.to('cpu')))

                test_loss += criterion(outputs, true_label).item()

        return true_label_list, pred_label_list, pred_score_list, test_loss

    def _add_roc_curve(self, fpr_list, tpr_list, tag):
        for class_idx in range(len(fpr_list)):
            fpr = fpr_list[class_idx] * 100
            tpr = tpr_list[class_idx] * 100

            for idx, x in enumerate(fpr):
                if np.isnan(x):
                    x = 0
                self.writer.add_scalars('ROC_AUC_curve/{}'.format(tag), {str(class_idx): tpr[idx]}, x)

    def _train(self, eps, optimizer, device, criterion, splits, batch_size, shuffle,
               log_per_n_epoch=5):
        # The exact implementation of train.
        logging.info('Action, Epoch,  Accuracy, Precision,    Recall,    ROCAUC,      Loss')
        train_idx, test_idx = splits
        trainloader = DataLoader(Subset(self.dataset, train_idx), batch_size=batch_size,
                                 shuffle=shuffle, num_workers=self.num_workers)
        testloader = DataLoader(Subset(self.dataset, test_idx), batch_size=batch_size,
                                shuffle=shuffle, num_workers=self.num_workers)

        n_batch = 1
        epoch_loss_list = []
        log_string = '{: >6}, {: >5n}, {: >9.2f}, {: >9.2f}, {: >9.2f}, {: >9.2f}, {: >9.4f}'
        tn_fpr, tn_tpr, tt_fpr, tt_tpr = None, None, None, None
        for epoch in range(eps):
            inputs = None
            running_loss = 0.0
            true_label_list, pred_label_list, pred_score_list = [], [], None

            for n_batch, data in enumerate(trainloader):
                inputs_type = torch.FloatTensor if device == 'cpu' else torch.cuda.FloatTensor

                # Get inputs and true label ready
                inputs = data[0].type(inputs_type).to(device)
                true_label = data[1][1].to(device)

                # Forward propagation
                outputs = self.net(inputs)

                # Calculate loss
                train_loss = criterion(outputs, true_label)

                # Backward propagation
                optimizer.zero_grad() # Reset parameter gradients
                # autograd.backward(train_loss, grad_tensor=None)
                train_loss.backward() # Backpropagate prediction loss
                optimizer.step()      # Adjust model parameters

                _, pred_label = torch.max(outputs.data, 1)
                pred_label_list.extend(pred_label.to('cpu'))
                true_label_list.extend(true_label.to('cpu'))

                if pred_score_list is None:
                    pred_score_list = outputs.data.to('cpu')
                else:
                    pred_score_list = torch.cat((pred_score_list, outputs.data.to('cpu')))

                running_loss += train_loss.item()

            running_loss /= n_batch

            if inputs is not None and epoch == 0:
                self.writer.add_graph(self.net, inputs)

            if (epoch + 1) % log_per_n_epoch == 0: # Log the evaluations per log_per_n_epoch.
                # Get training evaluations
                train_evm = self._eval_metric(true_label_list, pred_label_list, pred_score_list)
                tn_fpr, tn_tpr, tn_roc_auc, tn_prec, tn_accu, tn_rcll = train_evm

                # Write training evaluations to disk
                self.writer.add_scalar('ROC_AUC_score/Train', tn_roc_auc, epoch)
                self.writer.add_scalar('Precision/Train', tn_prec, epoch)
                self.writer.add_scalar('Accuracy/Train', tn_accu, epoch)
                self.writer.add_scalar('Recall/Train', tn_rcll, epoch)
                self.writer.add_scalar('Loss/Train', running_loss, epoch)
                logging.info(log_string.format("Train", epoch+1, tn_accu, tn_prec, tn_rcll,
                                               tn_roc_auc, running_loss))

                # Test
                test_output = self._test(self.net, testloader)
                true_label_list, pred_label_list, pred_score_list, tt_loss = test_output

                # Get test evaluations
                test_evm = self._eval_metric(true_label_list, pred_label_list, pred_score_list)
                tt_fpr, tt_tpr, tt_roc_auc, tt_prec, tt_accu, tt_rcll  = test_evm

                # Write test evaluations to disk
                self.writer.add_scalar('ROC_AUC_score/Test', tt_roc_auc, epoch)
                self.writer.add_scalar('Precision/Test', tt_prec, epoch)
                self.writer.add_scalar('Accuracy/Test', tt_accu, epoch)
                self.writer.add_scalar('Recall/Test', tt_rcll, epoch)
                self.writer.add_scalar('Loss/Test', tt_loss, epoch)
                logging.info(log_string.format('Test', epoch+1, tt_accu, tt_prec, tt_rcll,
                                               tt_roc_auc, tt_loss))

            epoch_loss_list.append(running_loss)

        self._add_roc_curve(tn_fpr, tn_tpr, 'Train')
        self._add_roc_curve(tt_fpr, tt_tpr, 'Test')
        self.train_loss_list.append(epoch_loss_list)

    def init(self):
        '''Init.'''
        # Check GPU or CPU
        self.writer = SummaryWriter(self.logging_path)
        self.net.to(self._check_device())
        return self

    def load_dataset(self, **kwargs):
        '''Load the dataset for train or test.'''
        logging.info("Start to load ASEDataset...")
        gene_id = kwargs.get('gene_id', self.gene_id)
        file_path_pool = kwargs.get('file_path_pool', self.file_path_pool)
        element_trans = kwargs.get('element_trans', SeqToHilbertAndMakeLabel())
        dataset_trans = kwargs.get('dataset_trans', MultipleTestAdjustment())

        if isinstance(gene_id, str):
            gene_id = [gene_id]
        elif not isinstance(gene_id, (tuple, list)):
            raise ValueError("gene_id should be a list/tuple/str!")

        # Chain sequences of mutliple gene ID.
        self.dataset = ASEDataset(file_path_pool=file_path_pool, gene_id=gene_id,
                                  element_trans=element_trans, dataset_trans=dataset_trans,
                                  num_workers=self.num_workers)
        logging.info("ASEDataset was ready...")

        return self

    def train_test_split(self, train_size=0.7):
        '''Split the dataset into train and test dataset.'''
        labels = list(self.dataset.get_labels())
        label_idx = range(len(labels))

        self.tt_splits = train_test_split(label_idx, train_size=train_size, stratify=labels)

        return self

    def train(self, eps=20, criterion=None, optimizer=None, learning_rate=1e-5, **kwargs):
        '''Train the model.

        Args:
            eps (int, optional, 20): Epoches.
            criterion (object, optional, None): Method to calculate loss.
            optimizer (object, optional, None): Method to optimize the model.
            learning_rate (float, optional, 1e-5): Learning rate. Works when `optimizer` is None.
            **kwargs (Any, optinal): Other keyword arguments.
        '''
        if criterion is None:
            criterion = nn.CrossEntropyLoss()

        if optimizer is None:
            optimizer = optim.Adam(self.net.parameters(), lr=learning_rate)

        batch_size = kwargs.get('batch_size', 32)
        shuffle = kwargs.get('shuffle', True)
        device = self._check_device()

        if self.tt_splits is None:
            logging.error('No train and test splits were found. Use train_test_split() first.')
            return self

        log_per_n_epoch = self.log_per_n_epoch
        logging.info("Start to train the model...")
        self._train(eps, optimizer, device, criterion, self.tt_splits, batch_size, shuffle,
                    log_per_n_epoch=log_per_n_epoch)
        logging.info("Model was ready!")

        self.writer.close()
        return self

    def loss_curve(self, loss_curve_path=None, **kwargs):
        '''Save the loss curve.'''
        train_loss_list = self.train_loss_list
        n_cv = len(train_loss_list)

        fig, axes = plt.subplots(n_cv)
        if n_cv == 1:
            loss_axe_pair = [[train_loss_list[0], axes]]
        elif n_cv > 1:
            loss_axe_pair = zip(train_loss_list, axes)
        else:
            logging.error('The loss list is empty!!')
            return self

        for loss, ax in loss_axe_pair:
            epoch_list = list(range(len(loss)))
            ax.plot(epoch_list, loss)
            ax.set_title('Loss per epoch')
            ax.set_xlabel('Epoch')
            ax.set_ylabel('Loss')

        fig.savefig(loss_curve_path, **kwargs)
        logging.info("Check the loss curve at: {}".format(loss_curve_path))

        return self

    def save_model(self, path, how='state'):
        '''Save the trained model, typlically the state dictionary.'''
        if self.net is None:
            logging.error('The model is empty, have you ever run the train() method?')
            sys.exit()

        if how == 'state':
            torch.save(self.net.state_dict(), path)
        elif how == 'model':
            torch.save(self.net, path)
        else:
            raise ValueError('Unsupported way to save the model: choose one from [state, model]')

        return self

    def test(self, model_state: str = None, testloader: DataLoader = None):
        '''Test for new dataset.'''
        self._test(model_state, testloader)

        return self


if __name__ == '__main__':
    logging.error('This file is meant to be imported.')
