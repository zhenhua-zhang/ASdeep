# Author : Zhenhua Zhang
# E-mail : zhenhua.zhang217@gmail.com
# Created: May 07, 2020
# Updated: Oct 06, 2021
"""A module to train a convolutionary neural network."""

import sys
import traceback
from argparse import Namespace

import h5py as h5
import numpy as np
from sklearn.model_selection import train_test_split as ttsplit
from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve

import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as func
from torch.utils.data import DataLoader
from torch.utils.data import Subset
from torch.utils.tensorboard import SummaryWriter

from .dataset import SubsetHilbertCurve
from .dataset import XyTransformer
from .dataset import ASEDataset
from .zutils import LogManager
from .zutils import pickup_model


class Trainer:
    """A class to train the neuroal network."""
    def __init__(self, model, dataset=None, log_output=None, log_n_epoch=5,
                 n_cpus=4, logman: LogManager = LogManager("Train")):
        self._model = model                 # NN model will be used
        self._n_cpus = n_cpus           # Number of workers to load data
        self._dataset = dataset         # Dataset to be loaded
        self._log_output = log_output   # Path for the TensorBoard logs
        self._log_n_epoch = log_n_epoch # Logging point per N epoches
        self._logman = logman           # Logging manager

        self._tt_splits = None     # train:test splits
        self._train_loss_list = [] # Loss

        self._writer = SummaryWriter(log_output)
        self._model.to(self.device)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, tb):
        if exc_type is not None:
            traceback.print_exception(exc_type, exc_value, tb)

        if hasattr(self._writer, "close"):
            self._writer.close()

    @property
    def model(self):
        return self._model

    @property
    def log_output(self):
        return self._log_output

    @property
    def device(self):
        return "cuda:0" if torch.cuda.is_available() else "cpu"

    @property
    def input_type(self):
        if self.device == "cpu":
            return torch.FloatTensor
        return torch.cuda.FloatTensor

    def _eval_matrix(self, y_true, y_pred, y_score, labels=[0, 1, 2],
                     scale=100):
        y_score = func.softmax(y_score, dim=1).data.numpy()
        fpr, tpr = [], []
        for idx, label in enumerate(labels):
            y_label_bin = [1 if x == label else -1 for x in y_true]
            _fpr, _tpr, _ = roc_curve(y_label_bin, y_score[:, idx])
            fpr.append(_fpr)
            tpr.append(_tpr)

        y_true_onehot = np.zeros((len(y_true), 3))
        for item_idx, class_idx in enumerate(y_true):
            y_true_onehot[item_idx][class_idx] = 1

        auc = 0
        try:
            auc = roc_auc_score(y_true_onehot, y_score, average="macro",
                                multi_class="ovr")
        except ValueError as e:
            self._logman.error(e)

        pre = precision_score(y_true, y_pred, average="macro", zero_division=1)
        rcl = recall_score(y_true, y_pred, average="macro", zero_division=1)
        acc = accuracy_score(y_true, y_pred)

        return fpr, tpr, auc * scale, pre * scale, rcl * scale, acc * scale

    def _test(self,
              model_state: str = None,
              testloader: DataLoader = None,
              criterion=None):
        if self._model is None:
            if model_state is not None:
                self._model.load_state_dict(torch.load(model_state))
            else:
                self._logman.error("The model state is missing")
                sys.exit()

        if testloader is None:
            self._logman.error("Missing testloader!! Exit...")
            return None, None, None, None

        if criterion is None:
            criterion = nn.CrossEntropyLoss()

        tt_y_true, tt_y_pred, tt_y_score = [], [], None
        self._model.to(self.device)
        test_loss = 0.0
        with torch.no_grad(): # The model should be trained.
            for data in testloader:
                matrix, label = data
                inputs = matrix.type(self.input_type).to(self.device) # inputs
                y_true = label.to(self.device) # true label

                outputs = self._model(inputs)

                _, y_pred = torch.max(outputs.data, 1)
                tt_y_pred.extend(y_pred.to("cpu"))
                tt_y_true.extend(y_true.to("cpu"))

                if tt_y_score is None:
                    tt_y_score = outputs.data.to("cpu")
                else:
                    tt_y_score = torch.cat((tt_y_score, outputs.data.to("cpu")))
                test_loss += criterion(outputs, y_true).item()

        return tt_y_true, tt_y_pred, tt_y_score, test_loss

    # Add ROC curve to the writer.
    def _add_roc_curve(self, fpr_list, tpr_list, tag):
        for class_idx in range(len(fpr_list)):
            fpr = fpr_list[class_idx] * 100
            tpr = tpr_list[class_idx] * 100

            for idx, x in enumerate(fpr):
                x = 0 if np.isnan(x) else x
                self._writer.add_scalars("ROC_AUC_curve/{}".format(tag),
                                         {str(class_idx): tpr[idx]}, x)

    # The exact implementation of train.
    def _train(self, n_epoch, optimer, criterion, splits, batch_size, shuffle,
               log_per_n_epoch=5):
        # import pdb; pdb.set_trace()
        hline = "{: >6},{: >5},{: >9},{: >9},{: >9},{: >9},{: >9}".format(
            "Act.", "Epo.", "Acc.", "Pre.", "Rec.", "ROCAUC", "Loss")
        self._logman.info("---- Training reports ----")
        self._logman.info(hline)
        train_idx, test_idx = splits
        trainloader = DataLoader(Subset(self._dataset, train_idx),
                                 batch_size=batch_size,
                                 shuffle=shuffle,
                                 num_workers=self._n_cpus)
        testloader = DataLoader(Subset(self._dataset, test_idx),
                                batch_size=batch_size,
                                shuffle=shuffle,
                                num_workers=self._n_cpus)

        epoch_loss_list = []
        fm = "{: >6},{: >5n},{: >9.2f},{: >9.2f},{: >9.2f},{: >9.2f},{: >9.4f}"
        tn_fpr, tn_tpr, tt_fpr, tt_tpr = None, None, None, None
        for epoch in range(1, n_epoch + 1):
            inputs = None
            running_loss = 0.0
            tn_y_true, tn_y_pred, tn_y_score = [], [], None

            for _, data in enumerate(trainloader):
                matrix, label = data
                inputs = matrix.type(self.input_type).to(self.device) # inputs
                y_true = label.to(self.device) # true label

                outputs = self._model(inputs) # Forward propagation
                train_loss = criterion(outputs, y_true) # Calculate loss

                optimer.zero_grad()   # Reset parameter gradients
                train_loss.backward() # Backpropagate prediction loss
                optimer.step()        # Adjust model parameters

                _, y_pred = torch.max(outputs.data, 1)
                tn_y_pred.extend(y_pred.to("cpu"))
                tn_y_true.extend(y_true.to("cpu"))

                if tn_y_score is None:
                    tn_y_score = outputs.data.to("cpu")
                else:
                    tn_y_score = torch.cat((tn_y_score, outputs.data.to("cpu")))

                running_loss += train_loss.item()

            if inputs is not None and epoch == 1:
                self._writer.add_graph(self._model, inputs)

            # Log the evaluations per log_per_n_epoch.
            if epoch % log_per_n_epoch == 0:
                # Get training evaluations
                tn_fpr, tn_tpr, tn_rauc, tn_pre, tn_rcl, tn_acc \
                        = self._eval_matrix(tn_y_true, tn_y_pred, tn_y_score)

                # Write training evaluations to disk
                self._writer.add_scalar("ROCAUCscore/Train", tn_rauc, epoch)
                self._writer.add_scalar("Precision/Train", tn_pre, epoch)
                self._writer.add_scalar("Accuracy/Train", tn_acc, epoch)
                self._writer.add_scalar("Recall/Train", tn_rcl, epoch)
                self._writer.add_scalar("Loss/Train", running_loss, epoch)
                self._logman.info(fm.format("Train", epoch, tn_acc, tn_pre,
                                            tn_rcl, tn_rauc, running_loss))

                # Test
                tt_y_true, tt_y_pred, tt_y_score, test_loss \
                        = self._test(self._model, testloader)

                # Get test evaluations
                tt_fpr, tt_tpr, tt_rauc, tt_pre, tt_rcl, tt_acc \
                        = self._eval_matrix(tt_y_true, tt_y_pred, tt_y_score)

                # Write test evaluations to disk
                self._writer.add_scalar("ROCAUCscore/Test", tt_rauc, epoch)
                self._writer.add_scalar("Precision/Test", tt_pre, epoch)
                self._writer.add_scalar("Accuracy/Test", tt_acc, epoch)
                self._writer.add_scalar("Recall/Test", tt_rcl, epoch)
                self._writer.add_scalar("Loss/Test", test_loss, epoch)
                self._logman.info(fm.format("Test", epoch, tt_acc, tt_pre,
                                            tt_rcl, tt_rauc, test_loss))
            epoch_loss_list.append(running_loss)

        self._add_roc_curve(tn_fpr, tn_tpr, "Train")
        self._add_roc_curve(tt_fpr, tt_tpr, "Test")
        self._train_loss_list.append(epoch_loss_list)

    def split_train_test(self, train_size=0.9):
        """Split the dataset into train and test dataset."""
        labels = list(self._dataset.get_labels())
        label_idx = range(len(labels))

        self._tt_splits = ttsplit(label_idx, train_size=train_size,
                                  stratify=labels)
        return self

    def train(self, epoches=20, criterion=None, optimizer=None,
              learning_rate=1e-5, batch_size=32, shuffle=True):
        """Train the model.

        Args:
            epoches: Epoches.
            criterion: Method to calculate loss.
            optimizer: Method to optimize the model.
            learning_rate: Learning rate. Works when `optimizer` is None.
            batch_size:
            shuffle:
        """
        if self._tt_splits is None:
            self._logman.error("No train/test splits. Use train_test_split().")
            return self

        if criterion is None:
            criterion = nn.CrossEntropyLoss()

        if optimizer is None:
            optimizer = optim.Adam(self._model.parameters(), lr=learning_rate)

        log_per_n_epoch = self._log_n_epoch
        self._logman.info("Start to train the model...")
        self._train(epoches, optimizer, criterion, self._tt_splits, batch_size,
                    shuffle, log_per_n_epoch=log_per_n_epoch)
        self._logman.info("Model was ready!")

        return self

    def save_model(self, path: str, how: str = "state", arch="model"):
        """Save the trained model, typlically the state dictionary."""
        if self._model is None:
            self._logman.error("Empty model, have you run the train method?")
            sys.exit()

        if how == "model":
            torch.save({arch: self._model}, path)
        else:
            torch.save({arch: self._model.state_dict()}, path)
            if how != "state":
                self._logman.warning(
                    "Unsupported way to save the model: choose one" +
                    " from [state, model]")

        self._logman.info("Check the model ({}) at {}".format(how, path))

        return self

    def test(self, model_state: str = None, testloader: DataLoader = None):
        """Test for new dataset."""
        self._test(model_state, testloader)

        return self


def train(args: Namespace, logman: LogManager = LogManager("Train")):
    """Train a CNN model."""
    database = args.database
    prebuilt_arch = args.prebuilt_arch
    model_state_path = args.model_state_path
    learning_rate = args.learning_rate
    epoches = args.epoches
    batch_size = args.batch_size
    n_base_pairs = args.n_base_pairs
    train_pp = args.train_pp
    log_output = args.log_output
    log_n_epoch = args.log_n_epoch
    random_state = args.random_state
    n_cpus = args.n_cpus

    np.random.seed(random_state)
    torch.manual_seed(random_state)
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True

    logman.info("---- Training parameters ----")
    logman.info(f"CPU counts:    {n_cpus}")
    logman.info(f"Batch size:    {batch_size}")
    logman.info(f"Random state:  {random_state}")
    logman.info(f"Epoch number:  {epoches}")
    logman.info(f"Architecture:  {prebuilt_arch}")
    logman.info(f"N base pairs:  {n_base_pairs}")
    logman.info(f"Training prop: {train_pp}")
    logman.info(f"Learning rate: {learning_rate}")

    trans = [SubsetHilbertCurve(n_base_pairs), XyTransformer()]
    net = pickup_model(prebuilt_arch)
    with h5.File(database, "r") as h5db:
        dataset = ASEDataset(database=h5db, transformers=trans, n_cpus=n_cpus)
        with Trainer(net, dataset, log_output, log_n_epoch, n_cpus) as trainer:
            (trainer.split_train_test(train_pp)
             .train(epoches, learning_rate=learning_rate,
                    batch_size=batch_size)
             .save_model(model_state_path, arch=prebuilt_arch))
