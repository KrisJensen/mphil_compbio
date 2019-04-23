#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
simulate mitochondrial inheritance using binomial distribution
"""
import matplotlib.pyplot as plt
import numpy as np

def run_sim(Ne = 3400, c=0.01):
    ns = [c]
    
    for i in range(2400):
        if Ne == 'exp':
            newNe = max([3400, 3400*np.exp(0.00629586756*(i-399))])
            n = np.random.binomial(newNe, ns[-1])/newNe
            
        else: n = np.random.binomial(Ne, ns[-1])/Ne
        ns.append(n)
        
    neander = False
    if ns[-1] > 0: neander = True
        
    return ns, neander

def plot_cgraph():
    
    N = 10000
    cs = np.arange(0.00, 1.02, 0.02)
    ps = [(1-np.sum([run_sim(Ne='exp', c=c)[1] for i in range(N)])/N) for c in cs]
    
    for i in range(len(cs)):
        print(cs[i], ps[i])

    plt.figure(figsize = (5,3))
    plt.plot(cs, ps, 'b-')
    plt.xlabel('c')
    plt.ylabel('$P(0)$')
    plt.savefig('figures/ccurve_exptemp.png', bbox_inches = 'tight', dpi=120)
    plt.show()
    plt.close()

def run_many_sims(nsim = 10000, Ne=3400, c = 0.01, fname='default'):
    
    if fname == 'default': fname = 'Ne'+str(Ne)+'_c'+str(c)[2:]
    
    n_neander = 0
    all_ns = np.zeros((nsim, 2401))
    tlost = []
    finalns = []
    for i in range(nsim):
        ns, neander = run_sim(Ne=Ne, c=c)
        if neander: n_neander+=1
        all_ns[i,:] = ns
        tlost.append(sum(np.array(ns) > 0)-1)
        finalns.append(ns[-1])
        
    ns = np.mean(all_ns, 0)
    
    plt.figure(figsize = (5,3))
    plt.plot(ns)
    plt.xlabel('generation')
    plt.ylabel('fraction neanderthal')
    plt.savefig('figures/'+fname+'_meann.png', bbox_inches = 'tight', dpi=120)
    plt.show()
    
    plt.figure(figsize = (5,3))
    plt.hist(tlost)
    plt.xlabel('time to loss / generations')
    plt.ylabel('frequency')
    plt.savefig('figures/'+fname+'_histloss.png', bbox_inches = 'tight', dpi=120)
    plt.show()
    
    bins = np.arange(0, 1.05, 0.05)
    plt.figure(figsize = (5,3))
    plt.hist(finalns, bins=bins)
    plt.xlabel('fraction Neanderthal at t=2400')
    plt.ylabel('frequency')
    plt.savefig('figures/'+fname+'_histfrac.png', bbox_inches = 'tight', dpi=120)
    plt.show()
    
    print(Ne, c, 'n neander', n_neander)
    
#run_many_sims(Ne=3400, c=0.01)
#run_many_sims(Ne=10000, c=0.01)
#run_many_sims(Ne='exp', c=0.01)

#run_many_sims(Ne=3400, c=0.10)
#run_many_sims(Ne=10000, c=0.10)
#run_many_sims(Ne='exp', c=0.1)

plot_cgraph()
    
        
    
    
    
    

