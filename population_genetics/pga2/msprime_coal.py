#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Code for simulating coalescents with msprime
"""

import numpy as np
import matplotlib.pyplot as plt
import msprime
import time
import pickle

nleaves = range(19, 38)
indps = [1.4833331078434333e-05,
0.00014833331078434333,
0.0007787498816178025,
0.002855416232598609,
0.008209321668721,
0.019702372004930403,
0.041046608343605004,
0.07622941549526645,
0.12863713864826212,
0.2001022156750744,
0.2901482127288579,
0.39565665372116987,
0.5110565110565111,
0.628992628992629,
0.7413127413127413,
0.8401544401544402,
0.918918918918919,
0.972972972972973,
1.0]
ps = []
for i in range(len(indps)):
    ps.append((1-sum(ps))*indps[i])

indps = [5.6167145457814216e-05,
0.0005055043091203279,
0.0024011454683215577,
0.008003818227738525,
0.02101002284781363,
0.04622205026518999,
0.08859226300828081,
0.15187245087133852,
0.23730070448646645,
0.3427676842582293,
0.46273637374860954,
0.5889372029527759,
0.7116324535679375,
0.8211143695014663,
0.9090909090909091,
0.9696969696969697,
1.0]
ps2 = []
for i in range(len(indps)):
    ps2.append((1-sum(ps2))*indps[i])

chr_ls=[249250621,
        243199373,
        198022430,
        191154276,
        180915260,
        171115067,
        159138663,
        146364022,
        141213431,
        135534747,
        135006516,
        133851895,
        115169878,
        107349540,
        102531392,
        90354753,
        81195210,
        78077248,
        63025520,
        59128983,
        51304566,
        48129895]
chr_ls.reverse()

Ne = 10**3
n = 38
nall = 31
tGen = 25

def get_nall(tree, ps=ps):
    #get most recent node with at least 31 descendants
    
    #talls = []
    #tmins = []
    
    nleaves = np.random.choice(range(len(ps),2*len(ps)), p=ps)
    
    for n in tree.nodes():
        nleaf = tree.get_num_leaves(n)
        if nleaf >= nleaves:
            return tree.time(n)
            #talls.append(tree.time(n))
        #elif nleaf == 2: tmins.append(tree.time(n))

    #return min(talls), min(tree.time(38))

def sim_genome(n=n,
               Ne=Ne,
               chr_ls = chr_ls,
               rho=0.5e-9*tGen, ps=ps):
    
    talls = []
    tmins = []
    for i, L in enumerate(chr_ls):
        print('\n starting new simulation')
        t = time.time()
        tree_sequence = msprime.simulate(sample_size=n,
                                 Ne=Ne,
                                 length=L,
                                 recombination_rate=rho)
        #seqs.append(tree_sequence)
        #trees = tree_sequence.aslist()
        print('simulation done', time.time()-t)
        
        tallst = []
        
        if n == 2:
            tallst.append(0)
            tmint = tree_sequence.get_time(n)
        else:
            for tree in tree_sequence.trees():
                tall = get_nall(tree, ps=ps)
                tallst.append(tall)
                tmint = tree_sequence.get_time(n)
        tallt = min(tallst)
        print('i='+str(i)+', L='+str(L)+', tall='+str(tallt)+
              ', tmin='+str(tmint)+', runtime='+str(time.time()-t))
        
        talls.append(tallt)
        tmins.append(tmint)
        
    #tall, tmin = min(talls), min(tmins)
    return talls, tmins

'''
def get_tpair(Ne=1000, nstudent = 17):
    tmins_all = []
    niter = int(nstudent*2*(nstudent-1)*2/2)
    niter = 50
    for i in range(niter): #number of unique pairs of chromosomes
        tmins = []
        t = time.time()
        tree_sequence = msprime.simulate(sample_size=2,
                                 Ne=Ne,
                                 length=2881033286,
                                 recombination_rate = 0.5e-9 * tGen)
        print('done simulating', i)
        tmin = tree_sequence.get_time(2)
        print('done trees', tmin, time.time()-t)
        tmins_all.append(tmin)
    return tmins_all
'''


def get_tpair(Ne=1000, nstudent = 17):
    tmins_all = []
    niter = int(nstudent*2*(nstudent-1)*2/2)
    niter = 50
    for i in range(niter): #number of unique pairs of chromosomes
        t = time.time()
        tree_sequence = msprime.simulate(sample_size=2,
                                 Ne=Ne,
                                 length=2881033286,
                                 recombination_rate = 0.5e-9 * tGen)
        print('done simulating', i)
        tmin = tree_sequence.get_time(2)
        print('done trees', tmin, time.time()-t)
        tmins_all.append(tmin)
    return tmins_all

def test_params_all():
    ts = []
    Nemin, Nemax = 500, 3000
    tgenmin, tgenmax = 20, 30
    rhomin, rhomax = 0.25e-9, 1e-9
    
    for i in range(100):
        tgen = np.random.uniform(tgenmin, tgenmax)
        Ne = np.random.uniform(Nemin, Nemax)
        rho = np.random.uniform(rhomin*tgen, rhomax*tgen)
        print('\n\n', tgen, Ne, rho)
        tall, tmin = sim_genome(Ne=Ne, rho=rho)
        ts.append(min(tall)*tgen)

    return ts

tree_sequence = msprime.simulate(sample_size=n, Ne=Ne)
tree = tree_sequence.first()
print(tree.draw(format="unicode"))
trees = tree_sequence.aslist()
print(get_nall(tree)) #get in generations

'''
ts = test_params_all()
pickle.dump(ts, open('ts_testparams.pickled2', 'wb'))
plt.figure(figsize = (5,3))
plt.hist(ts)
plt.xlabel('TMRCA / years')
plt.ylabel('frequency')
plt.savefig('figures/variableparams_allcoal2.png', dpi=120, bbox_inches = 'tight')
plt.show()
'''

tmins_all = []
for i in range(100):
    talls, tmins = sim_genome(Ne = 100000, n=34, ps=ps2, chr_ls = [10**7 for i in range(10)])
    tmins_all.append(tmins)

#newts = []
#for i in range(50):
#    print('\n\n\n new i', i)
#    tall, tmin = sim_genome()
#    newts.append(min(tall))

#plt.figure(figsize = (5,3))
#plt.hist([t*25 for t in newts])
#plt.xlabel('TMRCA / years')
#plt.ylabel('frequency')
#plt.savefig('figures/Ne1000_allcoal.png', dpi=120, bbox_inches = 'tight')

#pickle.dump([tall, tmin], open('fulltree_ts2.pickled', 'wb'))
#mins = get_tpair()

#newts = []
#for i in range(30):
#    print('\n\n\n new i', i)
#    tall, tmin = sim_genome(Ne = 10000, n=2)
#    newts.append(min(tmin))
#pickle.dump(newts, open('tmin_Ne10000.pickled', 'wb'))



#consider haploid genome
#tree_sequence = msprime.simulate(sample_size=n,
#                                 Ne=Ne,
#                                 length=100000000,
#                                 recombination_rate=0.5e-9*tGen)

#trees = tree_sequence.aslist()
    
#tree_sequence = msprime.simulate(sample_size=2,
#                                 Ne=1000,
#                                 length=115169878,
#                                 recombination_rate=0.5e-9*tGen)
#trees = tree_sequence.aslist()

