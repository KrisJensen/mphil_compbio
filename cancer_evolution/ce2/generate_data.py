#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Code for generating read data corresponding to the tree in assignemnt 1
"""

import numpy as np
import matplotlib.pyplot as plt
from numpy.random import negative_binomial, binomial
from scipy.stats import linregress
import random

#'Type' specifies which type of mutation to simulate
Type = 'after'
#Type = 'major'
#Type = 'minor'
random_type = True #pick a type at random
save = False #save figures

#parameters for negative binomial distribution
pbin = 0.3333
nbin = 50
#print(np.mean(negative_binomial(nbin, pbin, 10000)))
#print(np.var(negative_binomial(nbin, pbin, 10000)))
#print(negative_binomial(nbin, pbin, 10))

pi0 = 0.2 #healthy cell fraction
pis = [0.2, 0.2, 0, 0, 0.04, 0.12, 0.12, 0.12]
taus = [1, 0.75, 0.50, 0.30, 0.20, 0.15, 0.15, 0.15] #CCFs

nmuts = [5, 11, 3, 6, 6, 1, 5, 3]
varcounts, refcounts, fs, ccfs = [], [], [], [] #store data

Nsnv = 40
maj_cn = [1,2,3,4] # major copy numbers
maj_ps = [0.30, 0.30, 0.2,0.2] #probabilities
min_cn = [0,1,2] # minor copy numbers
min_ps = [1/4, 1/2, 1/4] # probabilities

#simulate copy number states
tmp1 = np.random.choice(maj_cn, Nsnv, p=maj_ps, replace = True)
tmp2 = np.random.choice(min_cn, Nsnv, p=min_ps, replace = True)
tots = tmp1 + tmp2
profile = [tmp1, tmp2, tots]


Write = True #write to file for future use
if Write:
    with open('PyClone_simdeep.tsv', 'w') as f:
                f.write('mutation_id\t'+
                        'ref_counts\t'+
                        'var_counts\t'+
                        'normal_cn\t'+
                        'minor_cn\t'+
                        'major_cn\n')
                
    with open('Ccube_simdeep.tsv', 'w') as f:
                f.write('mutation_id\t'+
                        'ccf_true\t'+
                        'minor_cn\t'+
                        'major_cn\t'+
                        'total_cn\t'+
                        'purity\t'+
                        'normal_cn\t'+
                        'mult_true\t'+
                        'vaf\t'+
                        'total_counts\t'+
                        'var_counts\t'+
                        'ref_counts\n')
    
nsite = -1
fifty = False#True
deep = True
if deep: nbin = 10*nbin

for i, tau in enumerate(taus):
    for j in range(nmuts[i]): #generate read counts for each site
        nsite += 1
        if random_type:
            if tmp2[nsite] == 0:
                ps = [0.5, 0.5, 0] #cannot observe what's not there
            else:
                ps = [0.5, 0.25, 0.25]
            Type = np.random.choice(['after', 'major', 'minor'], p = ps )
        if Type == 'after':
            nR = tots[nsite]
            nV = tots[nsite]
            nB = 1 #single allele
        elif Type == 'major':
            nR = 2
            nV = tots[nsite]
            nB = tmp1[nsite] #minor mutation
        elif Type == 'minor':
            nR = 2
            nV = tots[nsite]
            nB = tmp2[nsite] #major mutation
        else: print('something wrong')
        
        if fifty: #only seventy % CNVs
            if random.random() < 0.5:
                nR = 2
                nV = 2
                nB = 1
        
        #VAF
        p = (1-pi0)*tau*nB / (pi0*2 + (1-pi0)*(1-tau)*nR + (1-pi0)*tau*nV)
        r = negative_binomial(tots[nsite]*nbin,pbin) #total reads;
        R = binomial(r, p) #variant count
        print(R)
        refcount = r-R #reference allele
        varcounts.append(R) #store data
        refcounts.append(refcount)
        fs.append(p) #variant allele frequency
        ccfs.append(tau)
        
        if Write:
            with open('PyClone_simdeep.tsv', 'a') as f:
                f.write(str(nsite)+'\t'+
                        str(refcount)+'\t'+
                        str(R)+'\t'+
                        str(2)+'\t'+
                        str(tmp2[nsite])+'\t'+
                        str(tmp1[nsite])+'\n')
                
            with open('Ccube_simdeep.tsv', 'a') as f:
                f.write(str(nsite)+'\t'+
                        str(tau)+'\t'+
                        str(tmp2[nsite])+'\t'+
                        str(tmp1[nsite])+'\t'+
                        str(tots[nsite])+'\t'+
                        str(0.80)+'\t'+
                        str(2)+'\t'+
                        str(nB)+'\t'+
                        str(p)+'\t'+
                        str(r)+'\t'+
                        str(R)+'\t'+
                        str(refcount)+'\n')

#plot results
plt.figure(figsize = (4,3), dpi=120)
plt.bar(range(40), fs, align='edge')
plt.ylim(0,1.02)
plotmuts = np.cumsum([0]+nmuts)
for i in range(0, len(taus)):
    plt.plot([plotmuts[i], plotmuts[i+1]], [taus[i], taus[i]], 'k:')
plt.xlabel('locus')
plt.ylabel('VAF')
if save:
    plt.savefig(Type+'_bylocus.png')
plt.show()


plt.figure(figsize = (4,3), dpi=120)
s, i, r, p, std = linregress(ccfs, fs)
plt.plot([0, 1], np.array([0,1])*s+i, 'b:' )
plt.plot(ccfs, fs, 'bx', markersize=4)
plt.xlim(-0.05,1.05)
plt.ylim(-0.05,1.05)
plt.xlabel('CCF')
plt.ylabel('VAF')
plt.legend([str(np.round(s,2))+' CCF + '+str(np.round(i,2))+'   r='+str(np.round(r,2))])
if save:
    plt.savefig(Type+'_ccf_vaf.png')
plt.show()
