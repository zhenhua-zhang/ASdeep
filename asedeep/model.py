# -*- coding: utf-8 -*-
'''This is tiny CNN model copied from PyTorch tutorials.

Refer to: https://pytorch.org/tutorials/recipes/recipes/defining_a_neural_network.html
'''

import torch
import torch.nn as nn
import torch.nn.functional as fun

class ASEDeepNet(nn.Module):
    def __init__(self, input_size):
        super(ASEDeepNet, self).__init__()

        self.conv1 = nn.Conv2d(1, 32, 3, 1)
        self.conv2 = nn.Conv2d(32, 64, 3, 1)
        self.conv3 = nn.Conv2d(64, 128, 3, 1)
        self.conv4 = nn.Conv2d(128, 64, 3, 1)
        self.conv5 = nn.Conv2d(64, 16, 3, 1)
        self.conv6 = nn.Conv2d(16, 8, 3, 1)
        self.dropout1 = nn.Dropout2d(0.25)
        self.dropout2 = nn.Dropout2d(0.5)

        fc1_len = int(((input_size - 12) / 2)**2 * 8)
        self.fc1 = nn.Linear(fc1_len, 16)
        self.fc2 = nn.Linear(16, 3)

    def forward(self, x):
        x = self.conv1(x)
        x = fun.relu(x)

        x = self.conv2(x)
        x = fun.relu(x)

        x = self.conv3(x)
        x = fun.relu(x)

        x = self.conv4(x)
        x = fun.relu(x)

        x = self.conv5(x)
        x = fun.relu(x)
        
        x = self.conv6(x)
        x = fun.relu(x)

        x = fun.max_pool2d(x, 2)
        x = self.dropout1(x)
        x = torch.flatten(x, 1)

        x = self.fc1(x)
        x = fun.relu(x)
        x = self.dropout2(x)
        x = self.fc2(x)

        output = fun.log_softmax(x, dim=1)
        return output

