#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 12 13:22:44 2019

@author: kris
"""

import random
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

class net:
    def __init__(self, N=4, frac=0.5):
        
        self.G = nx.Graph()

        self.N = 0
        self.cols = {1:'b', -1:'r'}
        self.state = random.choices([-1,1], k=N)
        self.weights = np.zeros((N,N))
        
        #for i in range(N):
        #    self.add_node(random.sample([-1, 1], 1)[0])
        #    for j in range(i):
        #        self.add_con(str(i+1), str(j+1), random.sample([1,-1], 1)[0])

    def update_node(self, node):
        w = 0
        for j, props in self.G[node].items():
            #print(j.state, sign)
            w += props['weight']*self.G.nodes[j]['state']
            
        if w >= 0:
            self.G.nodes[node]['state'] = 1
            self.G.nodes[node]['colour'] = 'b'
        else:
            self.G.nodes[node]['state'] = -1
            self.G.nodes[node]['colour'] = 'r'
        
        return self.G.nodes[node]['state']
    
    def update_round(self):
        converged = True
        
        for node, props in self.G.nodes.items():
            if props['state'] != self.update_node(node):
                converged = False #have changed a node
            
        return converged
            
    def update_net(self):
        
        unconverged = True
        
        while unconverged:
            unconverged = not self.update_round()
            
    def add_node(self, state=1, cons = []):
        #self.nodes.append(node(state, cons))
        self.N+=1
        self.G.add_node(self.N, state = state, colour = self.cols[state])
        for j, sign in cons:
            self.G.add_edge(self.N, j, label = str(sign), sign = sign)
            
    def add_con(self, i, j, weight):
        
        self.G.add_edge(i, j, label=str(weight), weight=weight, colour = self.cols[weight])
        #self.cons[(node, j)] = weight
        #self.cons[(j, node)] = weight
            
    def add_cons(self, node, partners, weights):
        for i, j in enumerate(partners):
            if not (node, j) in self.cons.keys():
                weight = weights[i]
                self.add_con(node, j, weight)
                #col = self.cols[weight]
                #self.G.add_edge(node, j, label = str(weight), weight = weight, colour = col)
                #self.cons[(node, j)] = weight
                #self.cons[(j, node)] = weight
            
    def draw(self):
        ecols = [self.G[u][v]['colour'] for u,v in self.G.edges]
        ncols = [n['colour'] for n in self.G.nodes.values()]
        #print(cols)
        nx.draw(self.G, edge_color=ecols, node_color = ncols)
        plt.show()
        
    def drawsq(self):
        fig = plt.figure()
        ax = plt.subplots()[1]
        n = 0
        for i in range(int(self.N**0.5)):
            for j in range(int(self.N**0.5)):
                n+=1
                #state = self.G.nodes[str(n)]['state']
                rect = Rectangle((i, j), 1, 1, facecolor = self.G.nodes[n]['colour'])
                ax.add_patch(rect)
        for i in range(int(self.N**0.5)):
            plt.axvline(i, color='w')
            plt.axhline(i, color='w')
        plt.xlim([0, self.N**0.5])
        plt.ylim([0, self.N**0.5])
        plt.show()
        plt.close()
        
    def energy(self):
        E = 0
        for i in range(self.N):
            for j in range(i):
                n1 = i+1
                n2 = j+1
                E -= self.G[n1][n2]['weight'] *\
                    self.G.nodes[n1]['state'] *\
                    self.G.nodes[n2]['state']
        print(E)
        return E
    
    def learn_states(self, N, method = 'heb'):
        patterns = [np.array(random.choices([-1,1], k=self.N)) for n in range(N)]
        self.patterns = patterns
        self.learn_heb()
        
    def learn_heb(self):
        weights = np.zeros((self.N, self.N))
        for pattern in self.patterns:
            weights += pattern.reshape(self.N, 1).dot(pattern.reshape(1,self.N))
        for i in range(self.N): weights[i,i] = 0
        self.weights = weights
        
    def query(self, pattern, training = 'none'):
        for i in range(self.N):
            self.G.nodes[i]['state'] = pattern[i]
        self.update_net()
        newpat = np.array([self.G.nodes[i]['state'] for i in range(self.N)])
        
        if not training == 'none':
            error = 1 - newpat.dot(np.array(pattern))/self.N
            return error
        
        return [self.G.nodes[i]['state'] for i in range(self.N)]
            
    def __str__(self):
        out = ''
        for i, node in enumerate(self.nodes):
            out += (str(i)+' '+str(node.state)+'\n')
            
        return out
    
class node:
    def __init__(self, state=1, cons = []):
        self.state = state
        self.cons = []
        for j, sign in cons:
            self.add_con(j, sign)
        
        
    def add_con(self, j, sign = 1):
        self.cons.append((j, sign))
        
    def add_cons(self, js, signs):
        for i in range(len(js)):
            self.add_con(js[i], signs[i])
        

hop = net(N=16)


    
    
    
    

        
        
