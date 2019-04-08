#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 12:31:20 2019

@author: kris
"""

import numpy as np
import matplotlib.pyplot as plt
import copy

class perceptron:
    
    def __init__(self, N=10, gamma=-0.5):
        
        self.ws = np.zeros(N)
        self.gamma = gamma
        self.dim = N
        self.gamma = 0
        
        return
        
    
    def train(self, N = 100, method = 'heb', Print = False, task = 'sum'):
        
        train_data = np.random.choice([1,-1], (self.dim,N)) #training patterns
        if task == 'sum': vs = np.sign(np.sum(train_data,0)-self.gamma)
        elif task == 'product': vs = np.sign(np.prod(train_data, 0))
        else:
            print('task not recognized')
            return
        
        if Print: print(vs)
        self.data = train_data
        self.truth = vs
        
        if method == 'heb':
            self.train_heb()
        elif method == 'perceptron':
            self.train_perceptron()
        else:
            self.train_delta()

        if Print: print(self.ws)
        return
    
    def train_heb(self):
        self.ws = np.sum(self.data * self.truth, 1)
        return

    
    def train_perceptron_round(self, epsilon, thresh = 10**(-6), Print=False):
        oldws = copy.copy(self.ws)
        for i in range(self.data.shape[1]):
            u = self.data[:,i]
            self.ws += epsilon/2*(
                        self.truth[i] - 
                        np.sign(np.dot(self.ws,u)-self.gamma)) * u
            self.gamma -= epsilon/2*(
                        self.truth[i] - 
                        np.sign(np.dot(self.ws,u)-self.gamma))
                    
        err = np.sum(np.abs(self.ws-oldws))
        if Print: print('err:', err)
        if err < thresh: #converged
            return False
        else: return True #should still train
        
    def train_perceptron(self, epsilon = 0.01, nlim=1000):
        self.ws = np.zeros(self.dim)
        train = True
        n = 0
        while train and n < nlim:
            n += 1
            train = self.train_perceptron_round(epsilon)
        print('final gamma:', self.gamma)
        return
    
    def train_delta_round(self, epsilon, thresh = 10**(-6), Print=False):
        oldws = copy.copy(self.ws)
        delta = np.zeros(len(oldws))
        dgamma = 0
        for i in range(self.data.shape[1]):
            u = self.data[:,i]
            delta += (self.truth[i] - 
                        np.sign(np.dot(self.ws,u)-self.gamma)) * u
            dgamma += (self.truth[i] - 
                        np.sign(np.dot(self.ws,u)-self.gamma))
            
        self.ws += epsilon*delta
        self.gamma -= epsilon*dgamma
        err = np.sum(np.abs(self.ws-oldws))
        if Print: print('err:', err)
        if err < thresh: #converged
            return False
        else: return True #should still train
        
    def train_delta(self, epsilon = 0.01, nlim = 1000):
        self.ws = np.zeros(self.dim)
        train = True
        n = 0
        while train and n < nlim:
            n += 1
            train = self.train_delta_round(epsilon)
        print('final gamma:', self.gamma)
        return
        
    
    def query(self, u, Print = False):
        
        v = np.sign(np.dot(self.ws, u)-self.gamma)
        self.state = v
        if Print:
            print('sum:', np.sum(u))
            print('state:', v)
        return v
    
    def queries(self, nquery, plusses='random'):
        
        if plusses == 'random': queries = np.random.choice([-1, 1], (10,nquery))
        
        else:
            queries = np.zeros((10, nquery))
            pool = [1 for i in range(plusses)]+[-1 for i in range(10-plusses)]
            #print(pool)
            for i in range(nquery):
                queries[:,i] = np.random.choice(pool, 10, replace=False)
                #print(queries[:,i])
            
        
        truth = np.sign(np.sum(queries, 0)-self.gamma)
        vs = np.sign(np.dot(self.ws, queries)-self.gamma)
        performance = sum(np.array(vs) == np.array(truth))/nquery
        
        return performance
        
        '''
        positive = 0
        for i in range(nquery):
            u = np.random.choice([-1, 1], 10)
            truth = np.sign(np.sum(u)-self.gamma)
            v = self.query(u)
            if truth == v: positive += 1
        return positive/nquery
        '''
    
    def test_performance(self, method = 'heb', Ns = np.arange(1,400, 5), task = 'sum',
                         ntrial = 20, nquery=1000, Plot=True, plusses='random'):
        #plusses is nubmer of +1s in input
        performances = []
        for N in Ns:
            
            performance = []
            for j in range(10): #three separate nets
                self.train(N=N, method = method, task = task)
                for i in range(ntrial):
                    performance.append(self.queries(nquery, plusses=plusses))
            performance = np.mean(performance)
            performances.append(performance)
            
        if Plot:
            plt.figure()
            plt.plot(Ns, performances)
            plt.xlabel('N')
            plt.ylabel('performance')
            plt.ylim(0.4,1)
            plt.show()
            
        return Ns, performances
        
    def test_by_plus1(self, method = 'heb', fname='figures/test.png', task = 'sum'):
        plt.figure()
        nmax = 5
        cols = [[0,0,1], [0,1,0], [0,0,1], [0,1,1], [1,0,1], [1,1,0]]
        for nplus in range(nmax+1): #number of +1 vs -1
            print('new nplus:', nplus)
            Ns, performances = self.test_performance(method = method, Plot=False, plusses=nplus, task=task)
            plt.plot(Ns, performances, color = str(1-(nplus+1)/(nmax+1)), linestyle='-')
            print(sum(np.array(performances) < 0.95))
        plt.legend([str(val) for val in range(nmax+1)])
        
        plt.xlabel('N')
        plt.ylabel('performance')
        plt.ylim([0,1])
        plt.savefig(fname)
        plt.show()
        
    def compare(self, fname='figures/compare.png', task = 'sum'):
        
        plt.figure()
        cols = ['b', 'g', 'r']
        for i, method in enumerate(['heb', 'perceptron', 'delta']):
            print('new method:', method)
            Ns, performances = self.test_performance(method = method, Plot=False, task = task)
            plt.plot(Ns, performances, color = cols[i], linestyle='-')
            print(sum(np.array(performances) < 0.95))
        plt.ylim([0,1])
        plt.legend(['hebbian', 'perceptron', 'delta'])
        plt.xlabel('N')
        plt.ylabel('performance')
        plt.savefig(fname)
        plt.show()
        

p = perceptron()

p.compare(fname='figures/compare_sum.png')
p.compare(fname='figures/compare_prod.png', task='product')

#p.test_performance(method='delta')
p.test_by_plus1(method='delta', fname = 'figures/nplus_delta.png')


p = perceptron()
p.train(N=100, method='perceptron', task='product')        
p.query(np.random.choice([-1,1], 10), Print = True)
#p.test_performance(method='perceptron', task = 'product')
p.test_by_plus1(method='perceptron', fname = 'figures/nplus_perceptron.png')


#p = perceptron()
#p.train(N=100, method='heb')        
#p.query(np.random.choice([-1,1], 10), Print = True)
#p.queries(100, plusses=3)
#p.test_performance(method = 'heb')
p.test_by_plus1(method='heb', fname = 'figures/nplus_heb.png')



def plot_sum_prod():
    coords = [[1,1], [-1,-1], [1,-1], [-1,1]]
    
    plt.figure(figsize = (3,3))
    for coord in coords:
        if sum(coord) > 0:
            plt.plot(coord[0], coord[1], 'bo')
        elif sum(coord) < 0:
            plt.plot(coord[0], coord[1], 'ro')
        else:
            plt.plot(coord[0], coord[1], 'go')
            
    plt.xticks([-1,1], [-1,1])
    plt.yticks([-1,1], [-1,1])
    plt.xlim(-1.5,1.5)
    plt.ylim(-1.5,1.5)
    plt.savefig('figures/sum_2d.png', dpi = 120)
    plt.show()
    
    plt.figure(figsize = (3,3))
    for coord in coords:
        if coord[0]*coord[1] > 0:
            plt.plot(coord[0], coord[1], 'bo')
        elif coord[0]*coord[1] < 0:
            plt.plot(coord[0], coord[1], 'ro')
        else:
            plt.plot(coord[0], coord[1], 'go')
    plt.xticks([-1,1], [-1,1])
    plt.yticks([-1,1], [-1,1])
    plt.xlim(-1.5,1.5)
    plt.ylim(-1.5,1.5)
    plt.savefig('figures/prod_2d.png', dpi = 120)
    plt.show()
    
#plot_sum_prod()
    
    

    
    
    
    

