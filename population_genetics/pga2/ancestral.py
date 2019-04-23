#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Code for simulating coalescents without recombination
"""

import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

def add_sucs(G, node, newleaves):
    '''return list of leafnodes in order of their ancestors'''
    for newnode in G.successors(node):
        sucs = list(G.successors(newnode))
        if len(sucs) == 0: newleaves.append(newnode)
        else:
            newleaves = add_sucs(G, newnode, newleaves)
    return newleaves
    

def sort_graph(G):
    '''sort graph and annotate with initial x values to plot'''
    nodes = np.array(list(nx.topological_sort(G)))
    root = nodes[0]
    newleaves = add_sucs(G, root, [])
    xs = {}
    for i, leaf in enumerate(newleaves):
        xs[leaf] = i/(len(newleaves)-1)
    nodes = np.setdiff1d(nodes, np.array(newleaves))
    times= np.array([G.node[node]['time'] for node in nodes])
    args = np.argsort(times)
    nodes = nodes[args]
    times = times[args]
    nodelist = newleaves+list(nodes)
    
    times = [0 for i in range(len(nodelist)-len(times))]+list(times)
    #print(nodelist)
    return nodelist, xs, times

def get_xs(G, xs, nodelist):
    '''sort the nodes according to their ancestors so branches don't cross'''
    not_added = True
    while not_added:
        not_added = False
        for i, node in enumerate(nodelist):
            parents = list(G.predecessors(node))
            if len(parents) > 0:
                parent = parents[0]
                if not parent in xs.keys():
                    sucs = list(G.successors(parent))
                    for suc in sucs:
                        if not suc in xs.keys():
                            #print('dont have', suc)
                            xs = get_xs(G, xs, nodelist[i+1:])
                            not_added = True
                    #print(sucs)
                    xs[parent] = np.mean([xs[s] for s in sucs])
                        
    return xs


def draw_tree(G, relabel = False, fname = 'figures/drawtree_test.png'):
    '''plot the coalescent'''
    topo = list(nx.topological_sort(G))
    xs = {topo[0]: 0}
    nodelist, xs, times = sort_graph(G)
    xs = get_xs(G, xs, nodelist)
    plt.figure(figsize = (4,3))
    
    for node in nodelist:
        y1 = G.node[node]['time']
        x1 = xs[node]
        parents = list(G.predecessors(node))
        if len(parents) > 0:
            parent = parents[0]
            y2 = G.node[parent]['time']
            x2 = xs[parent]
            plt.plot([x1, x1], [y1, y2], 'k-')
            plt.plot([x1, x2], [y2, y2], 'k-')
        else:
            plt.plot([x1,x1], [y1, y1+0.1], 'k-')
        
    plt.xticks([], [])
    plt.ylim(0, max(times)+0.1)
    plt.xlim(-0.05, 1.05)
    plt.ylabel('time / N generations')
    if relabel: plt.savefig(fname, dpi = 120, bbox_inches = 'tight')
    else: plt.savefig(fname, dpi = 120, bbox_inches = 'tight')
    plt.show()
    

def sim_tree(n = 19, Print = True, dynamics = 'none'):
    '''simulate a coalescent with sample size n and N=1'''
    G= nx.DiGraph()
    G.add_nodes_from(range(1, n+1), time=0)
    
    t = 0 #start at time zero
    nnodes = n
    currnodes = [i for i in range(1,n+1)]
    
    while n > 1: #continue until everything has coalesced
        if dynamics == 'decreasing': #D > 0
            N = 5*(1/(1+np.exp(-2*(t-0.5)))-0.15)
        elif dynamics == 'increasing': #D < 0
            N = 10/(1+np.exp(2*t))
        else: N = 1 #D = 0
        
        rate = n*(n-1)/(2*N) #parameter of exponential distribution
        dt = np.random.exponential(scale = 1/rate) #next coalescence
        if Print: print(dt)
        t += dt #update time
        merge = np.random.choice(currnodes, 2, replace=False) #merge
        
        nnodes += 1
        G.add_node(nnodes, time=t) #update trees
        G.add_edge(nnodes, merge[0])
        G.add_edge(nnodes, merge[1])
        
        [currnodes.remove(node) for node in merge]
        currnodes.append(nnodes)
        n -= 1 #update number of lineages
        
    return G, t

def test_TMRCAs(ntrial = 100000, figname='figures/tmrca_hist.png'):
    #simulate ntrial trees and report TMRCA results
    tmrcas = []
    for i in range(ntrial):
        G, t = sim_tree(Print = False)
        tmrcas.append(t)
        
    plt.figure(figsize = (5,3))
    plt.hist(tmrcas, bins=50, density=True)
    plt.xlabel('TMRCA')
    plt.ylabel('frequency')
    plt.savefig(figname, dpi=120)
    print('mean:', np.mean(tmrcas), 'var:', np.var(tmrcas))
    
def plotdists():
    '''plot distributions of increasing and decreasing population sizes'''
    t = np.arange(0,3.1, 0.1)
    N1 = 10*(1/(1+np.exp(-2*(t-0.5))))-1.5
    N2 = 20/(1+np.exp(2*t))
    plt.plot(t, N1, 'b-')
    plt.plot(t,N2, 'g-')
    plt.xlabel('time / generations')
    plt.ylabel('N')
    plt.savefig('figures/changing_N.png')
    plt.show()
    
    
def plotsamps(n=19):
    '''plot the probability of sampling n individuals from K genomes'''
    Ks = range(n, 2*n+1)
    ps = [2**(2*n-K) * np.math.factorial(K) * np.math.factorial(n) / 
          (np.math.factorial(2*n) * np.math.factorial(K-n)) for K in Ks ]
    plt.figure(figsize = (5,3))
    plt.plot(Ks, ps, 'b-')
    plt.xticks( Ks, Ks )
    plt.xlabel('K')
    plt.ylabel('$p_{all}(K)$')
    plt.savefig('figures/testsamp_K.png', bbox_inches = 'tight', dpi=120)
    plt.show()
    for i in range(len(Ks)):
        print(Ks[i], ps[i])
    
plotsamps(19)
    
    
#plotdists()    

#G, t = sim_tree()
#draw_tree(G)
#test_TMRCAs()

#G, t = sim_tree(dynamics = 'decreasing')
#draw_tree(G, fname = 'figures/drawtree_decreasing2.png')

#G, t = sim_tree(dynamics = 'increasing')
#draw_tree(G, fname = 'figures/drawtree_increasing2.png')

#G, t = sim_tree(dynamics = 'constant')
#draw_tree(G, fname = 'figures/drawtree_constant2.png')

