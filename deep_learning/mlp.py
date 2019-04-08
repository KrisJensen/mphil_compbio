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
import copy
import time

class mlp():
    
    def __init__(self, N=20):
        #N is the size of the hidden layer
        #fully connected
        self.N = N
        
        data = pickle.load(open('MNIST/data.pickled', 'rb')) #load data
        train_ims, train_labs, test_ims, test_labs = data[0], data[4], data[1], data[5]
        
        self.ts = train_labs #store training data
        self.us = train_ims
        
        self.test_ts = test_labs #store test data
        self.test_us = test_ims
        
        ninput, nus = self.us.shape #number of inputs and input units
        
        self.ninput = ninput #store number of inputs
        self.ntest = self.test_us.shape[0] #store number of test datapoints
        
        self.w1s = np.random.normal(0, 1, (784, N)) #random initial weights
        self.w2s = np.random.normal(0, 1, (N, 10)) #random initial weights
        
        
    def train_round(self):
        '''trains a single epoch of the neural network using backpropagation'''
        for n in range(self.ninput): #for each image
            ts = self.ts[n,:] #target
            us = self.us[n,:] #input
            vs, zs = self.calc_vs_zs(us) #get output
            self.backpropagate(us, vs, zs, ts) #update weights
        
    def train(self, nepoch = 100, fname = 'test'):
        '''trains the network for nepoch epochs and plots a graph of performance'''

        n = 0
        performance = np.zeros(nepoch+1) #store performance        
        frac_corr0 = self.test_all()
        print('now at n='+str(n), 'performance='+str(np.round(frac_corr0,2)))
        performance[0] = frac_corr0
        
        
        while n < nepoch:#err > errlim and n < nmax:
            n += 1
            oldw1s = copy.copy(self.w1s)
            oldw2s = copy.copy(self.w2s)
            t = time.time()
            self.train_round()
            t = time.time() - t
            err = np.linalg.norm(self.w1s-oldw1s)+np.linalg.norm(self.w2s-oldw2s)
            frac_corr = self.test_all()
            print('now at n='+str(n), 'err='+str(np.round(err,2)),
                  'performance='+str(np.round(frac_corr,2)), 'time='+str(np.round(t,0)))
            performance[n] = frac_corr
            
        pickle.dump(performance, open(fname+'_performance.pickled', 'wb'))
        plt.figure(figsize = (4,2))
        plt.plot(range(nepoch+1), performance)
        plt.xlabel('epoch')
        plt.ylabel('performance')
        plt.xlim(0, nepoch)
        plt.ylim(0,1)
        plt.savefig('figures/'+fname+'_performance.png', dpi=120)
        plt.show()
            
        return performance
        
    def calc_vs_zs(self, us):
        '''for a given input, calculates hidden and output activities'''
        
        vs = np.tanh(np.dot(us, self.w1s)) #calculate hidden layer values
        zs = np.tanh(np.dot(vs, self.w2s)) #calculate output values
        
        return vs, zs
    
    def test_all(self):
        
        ts = np.argmax(self.test_ts, 1)
        correct = 0
        for i in range(self.ntest):
            target = ts[i]
            us = self.test_us[i,:]
            vs, zs = self.calc_vs_zs(us)
            prediction = np.argmax(zs) #use max function for prediction
            if prediction == target: correct += 1 #count number of correct predictions
        frac_corr = correct / self.ntest
        return frac_corr
        
    def backpropagate(self, us, vs, zs, ts, eta = 0.1):

        
        #update weights for second layer
        deltaks = np.zeros(10)
        #delta2s = np.zeros((self.N,10))
        dw2s = np.zeros((self.N, 10))
        for k in range(10): #for each output unit
            delta = (1 - zs[k]**2) * (ts[k]-zs[k]) #calculate delta
            deltaks[k] = delta
            for j in range(self.N): #for each hidden unit
                dw2s[j, k] = eta * delta * vs[j] #update rule
        self.w2s += dw2s
                
        #update weights for first layer
        deltajs = np.zeros(self.N)
        dw1s = np.zeros((784, self.N))
        for j in range(self.N): #for each hidden unit
            delta = (1-vs[j]**2) * np.sum(deltaks * self.w2s[j,:]) #calculate delta
            deltajs[j] = delta
            for i in range(784): #for each input unit
                dw1s[i,j] = eta * delta * us[i]
                
        self.w1s += dw1s

        return
        
m = mlp()

performance = m.train(fname='100_epoch', nepoch=100)
        

