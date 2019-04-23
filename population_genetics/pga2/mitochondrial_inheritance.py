#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Code for simulating mitochondrial inheritance

"""

import numpy as np
import matplotlib.pyplot as plt

def ang_brac(a, j):
    if j == 0: return 1
    else:
        P = 1
        for i in range(j):
            P *= (a-i)
        return P

def round_brac(a, j):
    if j == 0: return 1
    else:
        P = 1
        for i in range(j):
            P *= (a+i)
        return P

def g(n = 100, k=1, t=600, N=500):
    if k == 1:
        S = 1
        try:
            for i in range(2, n+1):
                rho = np.exp(-i*(i-1)* (t/N) / 2)
                ds = rho*(2*i-1)*((-1)**i) * ang_brac(n, i) / round_brac(n, i)
                S -= ds
                if ds == 0:
                    #print('change in sum is zero for g')
                    return S
                
                #print('new i:', i, 'S:', S)
        except OverflowError:
            print('overflowerror :( ')
            return S 
    else:
        S = 0
        try:
            for i in range(k, n+1):
                rho = np.exp(-i*(i-1)* (t/N) / 2)
                ds = (rho*(2*i-1)*((-1)**(i-k))*round_brac(k, i-1)*ang_brac(n,i)/ 
                    (np.math.factorial(k) * np.math.factorial(i-k) * round_brac(n, i)))
                S += ds
                
                if ds == 0:
                    #print('change in sum is zero for g')
                    return S
                #print('new i:', i, 'S:', S)
        except OverflowError:
            print('overflowerror :( ')
            return S  
    return S


def P0(c = 0.01, n=1000000, t=2400, N=10000):
    S = 0
    expected = 0
    sumgs = 0
    for k in range(1, n+1):
        newg = g(n=n, k=k, t=t, N=N)
        ds = ( newg * (1-c)**k )
        S += ds
        sumgs += newg
        expected += newg*k
        if ds == 0:
            #print('change in value is zero for P0')
            print('new exp g:', expected)
            #print('new sum g:', sumgs)
            return max(S, 0)
        #print('new k:', k, 'S:', S)
    return max(S, 0)

def plot_pgraph(n=986, t=100000, N=3400):
    t = t/20 #generations
    cs = np.arange(0.00, 1.01, 0.01)
    ps = [P0(c=c, n=n, t=t, N=N) for c in cs]
    plt.figure(figsize = (5,3))
    plt.plot(cs, ps, 'b-')
    plt.xlabel('c')
    plt.ylabel('$P(0)$')
    plt.savefig('figures/testpc_curve.png', bbox_inches = 'tight', dpi=120)
    plt.show()
    return cs, ps

def plot_Ngraph(n=1000000, t=60000, c=0.01):
    t = t/25 #generations
    Ns = np.arange(100, 10100, 100)
    ps = [P0(c=c, n=n, t=t, N=N) for N in Ns]
    plt.figure(figsize = (5,3))
    plt.plot(Ns, ps, 'b-')
    plt.xlabel('N')
    plt.ylabel('$P(0)$')
    plt.savefig('figures/popN_curve.png', bbox_inches = 'tight', dpi=120)
    plt.show()
    plt.close()
    return Ns, ps

def plot_nordborg():
    c1s, p1s = plot_pgraph(n=986, t=100000, N=3400)
    c2s, p2s = plot_pgraph(n=986, t=30000, N=3400)
    plt.figure(figsize = (5,3))
    plt.plot(c1s, p1s, 'b--')
    plt.plot(c2s, p2s, 'b-')
    plt.xlabel('c')
    plt.ylabel('$P(0)$')
    plt.legend(['$t_m=100000$', '$t_m=30000$'])
    plt.savefig('figures/testnordborg_curve.png', bbox_inches = 'tight', dpi=120)
    plt.show()
    plt.close()
    
def plot_cgraphs():
    c1s, p1s = plot_pgraph(n=1000000, t=60000, N=10000)
    c2s, p2s = plot_pgraph(n=1000000, t=60000, N=3400)
    plt.figure(figsize = (5,3))
    plt.plot(c2s, p2s, 'b-')
    plt.plot(c1s, p1s, 'b--')
    plt.xlabel('c')
    plt.ylabel('$P(0)$')
    plt.legend(['N=3400', 'N=10000'])
    plt.savefig('figures/testc_curve.png', bbox_inches = 'tight', dpi=120)
    plt.show()
    plt.close()
    
def plot_ngraphs(t=60000):
    t = t/25
    ns = range(8)
    ns = [10**n for n in ns]
    ps = []
    legends = []
    for N in [10000, 3400]:
        for c in [0.01, 0.1]:
            ps.append( [P0(c=c, n=n, t=t, N=N) for n in ns] )
            legends.append('N='+str(N)+', c='+str(c))
    
    plt.figure(figsize = (5,3))
    types = ['b--', 'y--', 'b-', 'y-']
    for i in range(len(ps)):
        plt.plot(ns, ps[i], types[i])
    
    plt.xscale('log')
    plt.xlabel('n')
    plt.ylabel('$P(0)$')
    plt.ylim(-0.01, 1.01)
    plt.legend(legends)
    plt.savefig('figures/sampn_curve.png', bbox_inches = 'tight', dpi=120)
    plt.show()
    plt.close()
    
    
#plot_ngraphs()
#plot_nordborg()
#plot_cgraphs()
#Ns, ps = plot_Ngraph()
#cs, ps = plot_pgraph()
    
'''
k1 = g()
k2 = g(k=2)
k3 = g(k=3)
k4 = g(k=4)

p1 = P0(n = 2)

p2 = P0(n = 100)

p3 = P0(n = 1000000)

p4 = P0(n = 1000000000)
'''





