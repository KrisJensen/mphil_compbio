#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 11:40:47 2019

@author: kris
"""

import matplotlib.pyplot as plt
import numpy as np
import networkx as nx


def topological_sort(G, relabel = False):
    nodes = [int(n) for n in G.nodes]
    inds = {}
    for i, n in enumerate(nodes): inds[n] = i
    
    done = False
    while not done:
        done = True
        for i in nodes:
            for j in nodes:
                if j in G.neighbors(i):
                    a = inds[i]
                    b = inds[j]
                    if a > b:
                        inds[i] = b
                        inds[j] = a
                        done = False
                        
    for node, ind in inds.items():
        nodes[ind] = node
        
    if not relabel: return nodes
    
    else:
        nodes.reverse()
        for i, n in enumerate(nodes):
            print(n, i)
            G.node[n]['newlabel'] = str(i+1)
            
        
        nodes.reverse()
        return nodes, G
    
    


def parse_tree(fname = "tree.dat"):
    G = nx.DiGraph()
    
    f = open(fname, 'r')
    for line in f:
        n, d1, l1, d2, l2 = line.split()
        G.add_node(int(n))
        G.add_edge(int(n), int(d1), length = float(l1))
        G.add_edge(int(n), int(d2), length = float(l2))
    f.close()
    
    #topo = topological_sort(G)
    topo, G = topological_sort(G, relabel = True)
    
    G.root = topo[0]
    G.node[G.root]['rootlength'] = 0
    G.node[G.root]['tier'] = 0
    
    
    for node in topo[1:]: #calculate length below root node; start from top
        pred = list(G.predecessors(node))[0]
        G.node[node]['rootlength'] = G[pred][node]['length'] + G.node[pred]['rootlength']
        G.node[node]['tier'] =  G.node[pred]['tier']+1
    
    topo.reverse() #reverse list
    
    for node in topo: #calculate number of nodes below each node; start from bottom
        ntot = 0
        nleaves = 0
        nint = 0
        length = 0
        for daughter in G.neighbors(node):
            ntot += G.node[daughter]['ntot']+1 #also add for this node
            nleaves += G.node[daughter]['nleaves']
            nint += G.node[daughter]['nint']
            if len(list(G.neighbors(daughter))) == 0: nleaves += 1 #daughter is a leaf
            else: nint += 1 #daughter is internal
            
        
        G.node[node]['ntot'] = ntot #total nodes below this node
        G.node[node]['nleaves'] = nleaves #leaves below this node
        G.node[node]['nint'] = nint #internal nodes below this node
    
    return G
    
def print_tree(G, relabel = False):
    
    nodes = np.array([n for n in G.nodes])
    ntots = np.array([G.node[n]['ntot'] for n in G.nodes])
    args = np.argsort(ntots)
    
    nodesn = nodes[args]
    ntots = ntots[args]
    #for n in nodesn:
    #    print(n,":", G.nodes[n]['ntot'], ' (', G.nodes[n]['nleaves'],
    #          ',', G.nodes[n]['nint'], ')')
        
    print('\n')
    
    if relabel:
        for n in nodesn:
            lab = G.node[n]['newlabel']
            print(lab,"&", n, '&', G.nodes[n]['ntot'], '&', G.nodes[n]['nleaves'],
                  '&', G.nodes[n]['nint'], '&', G.nodes[n]['rootlength'], '\\'+'\\')
        
    else:
        for n in nodesn:
            print(n,"&", G.nodes[n]['ntot'], '&', G.nodes[n]['nleaves'],
                  '&', G.nodes[n]['nint'], '&', G.nodes[n]['rootlength'], '\\'+'\\')
        
    print('\n')
        
    ls = np.array([G.node[n]['rootlength'] for n in G.nodes])
    args = np.argsort(-ls)
    nodesl = nodes[args]
    #for n in nodesl:
    #    print(n,"&", G.nodes[n]['rootlength'])



G = parse_tree()
#nx.draw(G, with_labels=True)    
print_tree(G)
print_tree(G, relabel = True)


#assume that ttwo changes do not occur over the course of a single step

def draw_tree(G, relabel = False):
    
    #topo = list(nx.topological_sort(G))
    topo = list(topological_sort(G))
    
    
    xs = {topo[0]: 0}
    
    plt.figure()
    if relabel:plt.text(-0.05, 0.08, str(G.node[G.root]['newlabel']))
    else: plt.text(-0.05, 0.08, str(G.root))
    
    slopes = [0.5, 0.4,0.3,0.2,0.1]
    
    for node in topo:
        y1 = -G.node[node]['rootlength']
        x1 = xs[node]
        daughters = list(G.neighbors(node))
        if len(daughters) == 1:
            y2 = -G.node[daughters[0]]['rootlength']
            x2 = x1
            xs[daughter[0]] = x2
            plt.plot([x1, x1], [y1, y2], 'b-')
            plt.text(x2+0.1, y2, str(daughters[0]))
        elif len(daughters) == 2:
            newxs = [xs[node]-1, xs[node]+1]
            newxs = [-1, 1]
            for i, daughter in enumerate(daughters):
                y2 = -G.node[daughter]['rootlength']
                slope = -1/(y2+1.5)
                x2 = x1+(y1-y2)*newxs[i]*slope#slopes[G.node[node]['tier']]
                #x2 = newxs[i]
                xs[daughter] = x2
                
                plt.plot([x1, x2], [y1, y2], 'b-')
                if relabel: s = G.node[daughter]['newlabel']
                else: s = str(daughter)
                if len(list(G.neighbors(daughter))) == 0:
                    plt.text(x2-0.05, y2-0.5, s)
                elif newxs[i] > 0:
                    plt.text(x2+0.03, y2, s)
                else:
                    plt.text(x2-0.16, y2, s)
        
    plt.xticks([], [])
    plt.ylim(-9, 0.5)
    if relabel: plt.savefig("figures/tree_relabelled.png")
    else: plt.savefig("figures/tree.png")
    plt.show()
    
draw_tree(G)
draw_tree(G, relabel = True)

    
    
    
