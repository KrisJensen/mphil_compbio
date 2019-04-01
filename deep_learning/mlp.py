#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 30 17:50:38 2019

@author: kris

Get mnist data from

https://www.python-course.eu/neural_network_mnist.php
"""

import numpy as np
import matplotlib.pyplot as plt
import pickle


data = pickle.load(open('MNIST/data.pickled', 'rb'))

train_ims, train_labs, test_ims, test_labs = data[0], data[4], data[1], data[5]


class mlp():
    
    def __init__(self, N=10):
        #N is the size of the hidden layer
        #fully connected
        self.N = N
        
        data = pickle.load(open('MNIST/data.pickled', 'rb'))
        train_ims, train_labs, test_ims, test_labs = data[0], data[4], data[1], data[5]
        
        self.ts = train_labs
        self.us = train_ims
        
        ninput, nus = self.us.shape #number of inputs and input units
        
        self.w1s = np.zeros((784, N))
        self.w2s = np.zeros((N, 10))
        
        
    def train_round(self):
        
        for n in range(self.N): #for each imge
            ts = self.ts[n,:] #target
            us = self.us[n,:] #input
            zs = calc_zs(us) #get output
            self.backpropagate(zs, ts) #update weights
        
    def train(self):
        while err < errlim:
            self.train_round()
        return
        
    def calc_zs(self, us):
        
        vs = np.tanh(np.dot(us, self.w1s)) #can I do this as single command?
        zs = np.tanh(np.dot(vs, self.w2s))
        
        return zs
        
    def backpropagate(self, zs, ts):
        
        
        delta2s = np.zeros((self.N,10))
        for k in range(10): #for each output unit
            for j in range(self.N): #for each hidden unit
                delta2s[j,k] = -(ts[k]-zs[k]) *  (1 - zs[k]**2)  * vs[j] #update rule
                
        for j in range(self.N): #for each hidden unit
            for i in range(784): #for each input unit
                #delta2s[j,k] = -(ts[k]-zs[k]) *  (1 - zs[k]**2)  * vs[j] #update rule
                print('ok')
        return
        
m = mlp()
        

