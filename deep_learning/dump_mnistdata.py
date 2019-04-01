#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 30 18:24:16 2019

@author: kris

dump data for MNIST
"""

import numpy as np
import matplotlib.pyplot as plt
import pickle

image_size = 28 # width and length
labels = 10 #  number of different labels
image_pixels = image_size ** 2
path = "MNIST/"
train_data = np.loadtxt(path + "mnist_train.csv", delimiter=",")
test_data = np.loadtxt(path + "mnist_test.csv", delimiter=",") 

train_imgs = np.asfarray(train_data[:, 1:]) / (255*0.99 + 0.01)
test_imgs = np.asfarray(test_data[:, 1:]) / (255*0.99 + 0.01)
train_labels = np.asfarray(train_data[:, :1])
test_labels = np.asfarray(test_data[:, :1])

lr = np.arange(labels)
# transform labels into one hot representation
train_labels_one_hot = (lr==train_labels).astype(np.float)
test_labels_one_hot = (lr==test_labels).astype(np.float)
# we don't want zeroes and ones in the labels neither:
train_labels_one_hot[train_labels_one_hot==0] = 0.01
train_labels_one_hot[train_labels_one_hot==1] = 0.99
test_labels_one_hot[test_labels_one_hot==0] = 0.01
test_labels_one_hot[test_labels_one_hot==1] = 0.99

with open("MNIST/data.pickled", "wb") as f:
    data = (train_imgs, 
            test_imgs, 
            train_labels,
            test_labels,
            train_labels_one_hot,
            test_labels_one_hot)
    pickle.dump(data, f)