#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Code for time-dependent selection questions
"""

import random
import numpy as np
import matplotlib.pyplot as plt
from numpy.random import poisson



def pfix(f, N=100):
    '''returns the probability of fixation given values of f and N'''
    return (1-np.exp(-2*f))/(1-np.exp(-2*N*f))


def plot_summary(ts, n1s, fixed, fname='test', xlim = [-0.05,0.1]):
    '''plot a series of graphs given values of n1 as a function of t'''
    
    #plot timeseries
    plt.figure()
    plt.plot(ts, n1s)
    plt.ylim(0,1000)
    plt.xlabel('time')
    plt.ylabel('$n_1$')
    plt.savefig('figures/'+fname+'_selection_trajec.png')
    plt.show()
    
    #plot histograms
    bins = np.linspace(min(fixed), max(fixed), 100)
    plt.figure()
    plt.hist(fixed, bins=bins)
    plt.plot([0,0], [0,200], 'k:', linewidth = 1)
    plt.xlim(xlim[0], xlim[1])
    plt.xlabel('$f$')
    plt.ylabel('frequency')
    plt.savefig('figures/'+fname+'_selection_hist2.png')
    plt.show()
    
    '''
    also plot predicted distributions
    ns = [n1, 1000-n1]
    xs = bins#np.linspace(-0.02, 0.08, 100)
    ps = [6*pfix(x) for x in xs]
    exps = [np.exp(-lambd*abs(x)) for x in xs]
    counts = [ ns[int((1+np.sign(x))/2)]/1000 for x in xs]
    
    combs = [20*exps[i]*ps[i]*counts[i] for i in range(len(xs))]
    
    plt.figure()
    plt.plot(xs, ps, 'r--')
    plt.plot(xs, exps, 'b--')
    plt.plot(xs, counts, 'c--')
    #plt.plot(xs, combs, 'g--')
    plt.legend(['$p_{fix}(f)$','$e^{-\lambda |f|}$','allele count'])
    plt.savefig('figures/'+fname+'_selection_dists.png')
    plt.show()
    
    plt.figure()
    plt.plot(xs, combs, 'g--')
    plt.savefig('figures/'+fname+'_combined_dists.png')
    plt.show()
    '''
    
def get_dists(xs, lambd=100):
    '''return theoretical distribution'''
    n1 = 820 #calculate analytically
    ns = [n1, 1000-n1]
    
    ps = np.array([pfix(x) for x in xs])
    ps /= np.amax(ps)
    
    exps = np.array([np.exp(-lambd*abs(x)) for x in xs])
    exps /= np.amax(exps)
    
    counts = np.array([ ns[int((1+np.sign(x))/2)]/1000 for x in xs])
    
    return ps, exps, counts
    
def plot_dists():
    '''plot expected n1 distributions in the different scenarios'''
    xs = np.linspace(-0.1, 0.1, 100)
    ps, exps, counts = get_dists(xs)
    
    plt.figure()
    plt.plot(xs, ps, 'r--')
    plt.plot(xs, exps, 'b--')
    plt.plot(xs, counts, 'c--')
    #plt.plot(xs, combs, 'g--')
    plt.ylim(0,1.05)
    plt.xlim(np.amin(xs), np.amax(xs))
    plt.legend(['$p_{fix}(f)$','$e^{-\lambda |f|}$','allele count'])
    plt.xlabel('$f$')
    plt.savefig('figures/selection_dists.png')
    plt.show()
    
    
    xs = np.linspace(-0.05, 0.1, 101)
    ps, exps, counts = get_dists(xs)
    combs = np.array([exps[i]*ps[i]*counts[i] for i in range(len(xs))])
    combs /= np.amax(combs)
    plt.figure()
    plt.plot(xs, combs, 'g--')
    plt.ylim(0,1.05)
    plt.xlim(np.amin(xs), np.amax(xs))
    plt.xlabel('$f$')
    plt.ylabel('frequency')
    plt.savefig('figures/all_combined_dists.png')
    plt.show()
    
    xs = np.linspace(-0.03, 0.1, 101)
    ps, exps, counts = get_dists(xs)
    combs_non = np.array([exps[i]*ps[i] for i in range(len(xs))])
    combs_non /= np.amax(combs_non)
    plt.figure()
    plt.plot(xs, combs_non, 'g--')
    plt.ylim(0,1.05)
    plt.xlim(np.amin(xs), np.amax(xs))
    plt.xlabel('$f$')
    plt.ylabel('frequency')
    plt.savefig('figures/no_n_combined_dists.png')
    plt.show()
    
    
    #with fixed fs, need to switch forward to switch back, so actual
    #switching interval is given by the mean switching intervals of the two times frequency
    #n1/n0 is given by the ratio of switching rates integrated over space
    xs = np.linspace(-0.035, 0.035, 100)
    L = len(xs)
    ps, exps, counts = get_dists(xs)
    sym = np.array([exps[i]/np.mean([1/ps[i], 1/ps[L-i-1]]) for i in range(len(xs))])
    sym /= np.amax(sym)
    plt.plot(xs, sym, 'g--')
    plt.ylim(0,1.05)
    plt.xlim(np.amin(xs), np.amax(xs))
    plt.xlabel('$f$')
    plt.ylabel('frequency')
    plt.savefig('figures/sym_combined_dists.png')
    plt.show()
    

    
    

def run_sim(N=100, lambd=100, changing = False, tlim=3000000,
            resample = False, fname='test', tau=5*10**(-6), xlim = [-0.05,0.1]):
    '''run simulation as described in the assinment.
    resample specifies which scenrio we consider, tau specifies the rate
    of change of the environment. tlim specifies the length of simulation'''
    
    change = 0 #start in default environment
    tswitch = poisson(1/tau) #timescale of switching
    
    f1s = np.random.exponential(1/lambd, 1000) #get an intial fitness at each locus
    f0s = -f1s #non-preferred allele
    
    states = np.ones(1000) #all loci start in state 1
    n1 = 1000 #we have 1000 loci in state 1
    
    t = 0 #initialize
    ts = [t]
    n1s = [n1]
    fixed = []
    
    while t < tlim:
        i = random.randint(0,999) #mutate random locus
        old = states[i]
        if old == 0: f = f1s[i] #switch to allele 1
        else: f = f0s[i] #switch to allele 0
        if resample: #draw new fitnesses
            if old == 0: f = np.random.exponential(1/lambd)
            else: f = -np.random.exponential(1/lambd)
            if change == 1: f = -f #change in environment

        p = pfix(f) #probability of fixation
        if random.random() < p: #fix allele
            new = 1-old
            states[i] = new
            n1 += new-old
            if t > 500000: #start collecting distribution after 500000 timesteps
                fixed.append(f)
                
        n1s.append(n1)
        dt = min(10, np.abs( 2*np.log(N-1)/(919*f) )) #fixation time
        t += dt
        ts.append(t)
        
        if changing:
            if t > tswitch:
                #change environment
                f1s = -f1s
                f0s = -f0s
                change = 1-change
                tswitch += poisson(1/tau) #update switching time
        
    plot_summary(ts, n1s, fixed, fname=fname, xlim = xlim)
    
        
#plot_dists() #plot expected distributions

#run simulations with constant environment
#run_sim(fname='constant', xlim = [-0.035, 0.035]) #scenario 2
#run_sim(resample = True, fname='resample', xlim=[-0.05,0.1]) #scenario 1

#run simulations with changing environments
#run_sim(fname='changing', changing=True, xlim = [-0.03, 0.1])
#run_sim(fname='changing_resample', changing=True, resample = True, xlim=[-0.03,0.1])

    

#run_sim(fname='tau5e7', changing=True, tau=5*10**(-7), resample=False, xlim = [-0.03, 0.1])
#run_sim(fname='tau5e5', changing=True, tau=5*10**(-6.5), resample=False, xlim = [-0.03, 0.1])
#run_sim(fname='tau5e3', changing=True, tau=5*10**(-3), resample=False, xlim = [-0.03, 0.1])  

#questions: is the probability of fixation per unit time or absolute?
#how do we have a negative time to fixation?
    
    


