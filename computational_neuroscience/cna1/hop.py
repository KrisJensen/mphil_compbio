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
import copy
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
import pickle
import scipy.stats

random.seed(22021039)

class net:
    def __init__(self, N=4, frac=0.5):
        
        self.G = nx.Graph()
        self.N = N
        self.cols = {1:'b', -1:'r'}
        self.state = np.array(random.choices([-1,1], k=N))
        self.weights = np.zeros((N,N))


    def update_node(self, node):
        w = 0
        w = self.state.dot(self.weights[node, :])
            
        if w >= 0:
            self.state[node] = 1
        else:
            self.state[node] = -1
        
        return self.state[node]
    
    def update_round(self):
        converged = True
        for node in range(self.N):
            if self.state[node] != self.update_node(node):
                converged = False #have changed a node
        return converged
            
    def update_net(self):
        
        unconverged = True
        while unconverged:
            unconverged = not self.update_round()

    def draw(self):
        plt.figure()
        ax = plt.subplots()[1]
        n = 0
        for i in range(int(self.N**0.5)):
            for j in range(int(self.N**0.5)):
                n+=1
                #state = self.G.nodes[str(n)]['state']
                rect = Rectangle((i, j), 1, 1, facecolor = self.cols[self.state[n]])
                ax.add_patch(rect)
        for i in range(int(self.N**0.5)):
            plt.axvline(i, color='w')
            plt.axhline(i, color='w')
        plt.xlim([0, self.N**0.5])
        plt.ylim([0, self.N**0.5])
        plt.show()
        plt.close()
        
    def energy(self):
        E = -0.5*self.state.dot(self.weights).dot(self.state)
        print(E)
        return E
    
    def string_to_binary(self, st):
        nmiss = int(self.N/7-len(st))
        #print(nmiss)
        binary = []
        for x in st:
            if x == ' ': binary += [1, 0, 1, 0, 1, 0, 1]
            else: binary += [int(char) for char in format(ord(x), 'b')]
        for n in range(nmiss):
            binary += [[0,1,0,1,0,1,0], [1,0,1,0,1,0,1]][np.random.choice([0,1])]
        binary = np.array(binary)
        binary[binary == 0] = -1
        return(binary)
        
    def binary_to_string(self, binary):
        binary[binary == -1] = 0
        while (np.array_equal(binary[-7:], [0,1,0,1,0,1,0]) or
               np.array_equal(binary[-7:], [1,0,1,0,1,0,1])):
            binary = binary[:-7]
            
        binary = [str(Bin) for Bin in binary]
        byte = []
        for i in range(int(len(binary)/7)):
            newbin = ''.join(binary[:7])
            #print(newbin)
            if newbin == '1010101': byte.append(32)
            else: byte.append(int(newbin, 2))
            #print(byte)
            binary = binary[7:]
            #print(binary)
            
        text = bytes(byte).decode('utf-8')
        
        return(text)
        
    
    def learn_strings(self, strings = ['hello world']):
        patterns = []
        for st in strings:
            patterns.append(self.string_to_binary(st))
            
        self.patterns = patterns
        
        self.learn_opt()
        
    def query_string(self, st):
        pattern = self.string_to_binary(st)
        self.state = pattern
        self.update_net()
        newst = self.binary_to_string(self.state)
        print(newst)
        return newst
        
    
    def learn_states(self, N, method = 'heb', alpha=0.5, floss=0):
        #alpha is proportion of states that are +1
        #patterns = [np.array(random.choices([-1,1], k=self.N)) for n in range(N)]
        patterns = [np.array(np.random.choice([-1,1],
                    size=self.N, p=[1-alpha, alpha])) for n in range(N)]
        self.patterns = patterns
        if method == 'heb': self.learn_heb(floss=floss)
        elif method == 'opt': self.learn_opt(floss=floss)
        
    def learn_heb(self, floss=0):
        weights = np.zeros((self.N, self.N))
        data = np.array(self.patterns)
        self.data=data
        weights = data.transpose().dot(data)
        #for pattern in self.patterns:
        #    weights += pattern.reshape(self.N, 1).dot(pattern.reshape(1,self.N))
        for i in range(self.N): weights[i,i] = 0
        
        if floss > 0:
            indices = np.random.choice(range(int((self.N**2-self.N)/2)), replace=False,
                           size=int((self.N**2-self.N)/2 * floss))
            inds = []
            for i in range(self.N):
                for j in range(i):
                    inds.append([i, j])
            for ind in indices:
                weights[inds[ind][0], inds[ind][1]] = 0
                weights[inds[ind][1], inds[ind][0]] = 0
                
        self.weights = weights
        
    def learn_opt(self, floss=0, nmax=100, eta = 0.1):
        self.learn_heb(floss = floss) #initialize with hebbian weights
        t = copy.deepcopy(self.data)
        t[t==-1] = 0 #need for errors
        derr, err0 = 1, np.inf
        n=0
        while n < nmax and derr > 0:
            n+=1
            a = self.data.dot(self.weights)
            e = self.data-np.sign(a)
            err = np.mean(np.abs(e))
            #print('err:', err)
            derr = err0-err
            err0 = err
            gw = self.data.transpose().dot(e)
            gw = gw + gw.transpose() #symmetrize
            self.weights = self.weights + eta*gw
            #print('ok')
        
    def query(self, pattern, training = [], alpha=0.5):
        self.state = np.array(pattern)
        self.update_net()        
        if len(training) == self.N:
            s = alpha**2+(1-alpha)**2-2*alpha*(1-alpha) #error term from sparseness
            dp = np.abs(self.state.dot(np.array(training))) #negative patterns are also allowed
            error = 1/(1-s) - dp/self.N / (1-s)
            #print('error:', error)
            return error
        
        return self.state
    
    def test_states(self, nnets=5, npats = 1000, ferrs=0.1, nsamps=100,
                    alpha=0.5, floss=0, method='heb'):
        nerrs = int(ferrs*self.N)
        meanerrs = []
        for i in range(nnets):
            errs = []
            self.learn_states(npats, alpha=alpha, floss=floss, method=method)
            for j in random.choices(range(npats)):
                newpat =copy.deepcopy(self.patterns[j])
                ks = random.sample(range(self.N),nerrs)
                for k in ks: newpat[k] *= (-1)
                err = self.query(newpat, training=self.patterns[j], alpha=alpha)
                errs.append(err)
            meanerr = np.mean(errs)
            #print('mean error:', meanerr)
            meanerrs.append(meanerr)
            
        #print(meanerrs)
        return(np.mean(meanerrs))
            
    def __str__(self):
        out = ''
        for i, node in enumerate(self.nodes):
            out += (str(i)+' '+str(node)+'\n')
            
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
            
            
def scan_vars(Ns = [10, 600, 10], npats = [1, 200, 5], fname = 'scan_N_npat.png', method='heb'):
    
    Ns = np.arange( Ns[0], Ns[1], Ns[2] )
    npats = np.arange( npats[0], npats[1], npats[2] )
    vals = np.zeros((len(Ns), len(npats)))
    for i, N in enumerate(Ns):
        for j, npat in enumerate(npats):
            hop = net(N=int(N))
            err = hop.test_states(nnets=10, npats=int(npat), ferrs=0.1,
                                  nsamps=100, alpha=0.5, method=method)
            print(N, npat, err)
            vals[i,j] = err
    thresh = 0.138*Ns
            
    plot_3d(Ns, npats, vals, xlab='N', ylab='patterns',
            zlab='$\epsilon$',fname=fname, vmin=-0.4, addline=[Ns, thresh, 0.35])
            
    return Ns, npats, vals

def scan_alphas(alphas = [0.5,1.0, 0.01], npats = [1, 100, 2], N = 300,
                fname = 'test.png', method = 'heb'):
    
    alphas = np.arange( alphas[0], alphas[1], alphas[2] )
    npats = np.arange( npats[0], npats[1], npats[2] )
    vals = np.zeros((len(alphas), len(npats)))
    for i, alpha in enumerate(alphas):
        for j, npat in enumerate(npats):
            hop = net(N=int(N))
            err = hop.test_states(nnets=10, npats=int(npat), ferrs=0.1,
                                  nsamps=100, alpha=alpha, method=method)
            print(alpha, npat, err)
            vals[i,j] = err
    plot_3d(alphas, npats, vals, xlab="alpha", ylab='patterns',
            zlab="$\epsilon$",fname='scan_alpha_npat.png', vmin=-0.4, view = [55,250])
            
    return alphas, npats, vals

def scan_ferrs(ferrs = [0.0,0.5, 0.01], npats = [1, 100, 2], N = 300,
               fname = 'test.png', method='heb'):
    
    ferrs = np.arange( ferrs[0], ferrs[1], ferrs[2] )
    npats = np.arange( npats[0], npats[1], npats[2] )
    vals = np.zeros((len(ferrs), len(npats)))
    for i, ferr in enumerate(ferrs):
        for j, npat in enumerate(npats):
            hop = net(N=int(N))
            err = hop.test_states(nnets=10, npats=int(npat), ferrs=ferr,
                                  nsamps=100, alpha=0.5, method=method)
            print(ferr, npat, err)
            vals[i,j] = err
    plot_3d(ferrs, npats, vals, xlab="$f_{err}$", ylab='patterns',
            zlab="$\epsilon$",fname='scan_ferr_npat.png', vmin=-0.4, view = [55,250])
            
    return ferrs, npats, vals

def scan_floss(floss = [0.0,0.99, 0.02], npats = [1, 100, 2], N = 300, fname = 'test.png'):
    
    floss = np.arange( floss[0], floss[1], floss[2] )
    npats = np.arange( npats[0], npats[1], npats[2] )
    vals = np.zeros((len(floss), len(npats)))
    for i, flo in enumerate(floss):
        for j, npat in enumerate(npats):
            hop = net(N=int(N))
            err = hop.test_states(nnets=10, npats=int(npat), ferrs=0.1, nsamps=100,
                                  alpha=0.5, floss=flo)
            print(flo, npat, err)
            vals[i,j] = err
            
    plot_3d(floss, npats, vals, xlab="$f_{loss}$", ylab='patterns',
            zlab="$\epsilon$",fname='scan_floss_npat.png', vmin=-0.4,
            view = [55,250], alpha=0.7)
            
    return ferrs, npats, vals

def plot_3d(xs, ys, vals, xlab='', ylab='', zlab='', fname='test.png',
            vmin=0, addline=[], view = [55,310], alpha=0.6, zlim = []):
    xsa = np.repeat(xs,len(ys)).reshape(len(xs),len(ys))
    ysa = np.transpose(np.repeat(ys,len(xs)).reshape(len(ys),len(xs)))
    
    fig = plt.figure(figsize = (7,4))
    ax = fig.add_subplot(111, projection='3d')
    
    print(xsa.shape, ysa.shape, vals.shape)
    
    if len(addline) > 0:
        ax.plot(addline[0], addline[1], addline[2], 'k-', alpha=1, zorder=2)
    
    surf = ax.plot_surface(xsa, ysa, vals, cmap=cm.Greys, antialiased = False,
                           vmin = vmin, vmax = np.amax(vals), zorder=1, alpha=alpha)
    
        
    fig.colorbar(surf, shrink=0.5, aspect=6, pad=0.05)
    
    #make things look nice
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_zlabel(zlab)
    plt.xlim(xs[0], xs[-1])
    plt.ylim(ys[0], ys[-1])
    if len(zlim)>0: ax.set_zlim(zlim[0], zlim[1])
    #save and show
    ax.view_init(view[0], view[1])#55, 310)
    ax.dist = 10.5
    
    plt.savefig(fname, dpi=360)
    plt.show()
    
def test_learns(Ns = [10, 600, 10], npats = [1, 200, 5], fname='compare.png'):
    
    N, npat, vals1 = scan_vars(
        Ns = Ns, npats = npats,
              fname = 'test.png', method='heb')
    
    N, npat, vals2 = scan_vars(
        Ns = Ns, npats = npats,
              fname = 'test.png', method='opt')
    vals = vals1-vals2
    plot_3d(N, npat, vals, xlab="N", ylab='patterns',
            zlab="$\epsilon_{heb} - \epsilon_{opt}$",fname='scan_N_npat_hebopt.png', vmin=0,
            view = [55,250], alpha=1)
    
    return N, npat, vals1, vals2

def get_thresh(vals, Ns, nsamps, inf = True):
    infsamps = []
    
    if inf:
        for i in range(vals.shape[0]):
            maxdif = 0
            olddif = 0
            for j in range(1, vals.shape[1]):
                newdif = vals[i,j]-vals[i,j-1]
                if newdif > maxdif and olddif >= 0:
                    maxdif = newdif
                    maxsamp = (nsamps[j]+nsamps[j-1])/2
                olddif=newdif
            infsamps.append(maxsamp)
        
        
    else:
        for i in range(vals.shape[0]):
            maxdif = 0
            for j, val in enumerate(vals[i,:]):
                if val > 0.2:
                    if j == 0: infsamps.append(0)
                    else: infsamps.append((nsamps[j]+nsamps[j-1])/2)
                    break
                
    s, i, r, p, std = scipy.stats.linregress(Ns, infsamps)
        
    return Ns, infsamps, s, i, r

def plot_infls(vals1, Ns, nsamps, fname='compare', vals2=[], xlab = 'N', line=True):
    
    names = ['inflection', 'thresh']
    for i, inf in enumerate([True, False]):
        N, inf1, s1, i1, r1 = get_thresh(vals1, Ns, nsamps, inf = inf)
        if len(vals2) > 0: N, inf2, s2, i1, r1 = get_thresh(vals2, Ns, nsamps, inf=inf)
        plt.figure()
        plt.plot(Ns, inf1, 'b--')
        if len(vals2)>0:
            plt.plot(Ns, inf2, 'r--')
            plt.legend(['Hebbian (a='+str(np.round(s1,3))+')',
                        'Optimal (a='+str(np.round(s2,3))+')'])
            plt.plot(Ns, 0.138*Ns, 'k:')
        else:
            if line:
                plt.plot([N[0], Ns[-1]], [Ns[0]*s1+i1, Ns[-1]*s1+i1], 'b:')
                plt.plot(Ns, 0.138*Ns, 'k:')
                plt.legend(['Data',
                        'Fit (a='+str(np.round(s1,3))+' r='+str(np.round(r1,3))+')',
                        '0.138N'])

        plt.xlabel(xlab)
        plt.ylabel('patterns')
        plt.savefig(fname+'_'+names[i]+'.png')
        plt.show()
        
def test_err_loss(floss = [0,0.99,0.02], ferrs = [0,0.5,0.01], npats = [1,int(0.25*300),4]):
    N = 300
    floss = np.arange( floss[0], floss[1], floss[2] )
    ferrs = np.arange( ferrs[0], ferrs[1], ferrs[2] )
    npats = np.arange( npats[0], npats[1], npats[2] )
    vals = np.zeros((len(floss), len(ferrs)))
    for i, flo in enumerate(floss):
        tempvals = np.zeros((len(ferrs), len(npats)))
        for j, ferr in enumerate(ferrs):
            for k, npat in enumerate(npats): 
                hop = net(N=int(N))
                err = hop.test_states(nnets=10, npats=npat, ferrs=ferr, nsamps=100,
                                  alpha=0.5, floss=flo)
                print(flo, ferr, npat, err)
                tempvals[j,k] = err
                
        #print(tempvals)
        #print(tempvals.shape, ferrs.shape, npats.shape)
        Ns, infsamps, s, intercept, r = get_thresh(tempvals, ferrs, npats, inf = False)
        vals[i, :] = np.array(infsamps)/300
        
    plot_3d(floss, ferrs, vals, xlab="$f_{loss}$", ylab='$f_{err}$',
            zlab="capacity",fname='test_err_loss.png', vmin=0,
            view = [55,250], alpha=1)
        
        
    return floss, ferrs, vals
        
    
'''        
    

#hop = net(N=100)
#hop.test_states(nnets=3, npats = 10, ferrs=0.1, nsamps=100, alpha=0.1)

#st = 'test'
#print(' '.join(format(ord(x), 'b') for x in st))

#Ns, npats, vals = scan_vars()
#pickle.dump([Ns, npats, vals], open('scan_Ns_npats.pickled', 'wb'))
Ns, npats, vals = pickle.load(open('scan_Ns_npats.pickled', 'rb'))
thresh = 0.138*Ns
plot_3d(Ns, npats, vals, xlab='N', ylab='patterns',
            zlab="$\epsilon$",fname='scan_N_npat.png', vmin=0, alpha=1)
            #addline = [Ns,thresh,0.35] )
plot_infls(vals, Ns, npats, fname='hebb')



#alphas, npatsa, valsa = scan_alphas()
#pickle.dump([alphas, npatsa, valsa], open('scan_alphas_npats.pickled', 'wb'))
alphas, npatsa, valsa = pickle.load(open('scan_alphas_npats.pickled', 'rb'))
plot_3d(alphas, npatsa, valsa, xlab="alpha", ylab='patterns',
            zlab="$\epsilon$",fname='scan_alpha_npat.png',
           vmin=0, view = [55,250], alpha=1)

#ferrs, npatsf, valsf = scan_ferrs()
#pickle.dump([ferrs, npatsf, valsf], open('scan_ferrs_npats.pickled', 'wb'))
ferrs, npatsf, valsf = pickle.load(open('scan_ferrs_npats.pickled', 'rb'))  
plot_3d(ferrs, npatsf, valsf, xlab="$f_{err}$", ylab='patterns',
            zlab="$\epsilon$",fname='scan_ferr_npat.png', vmin=0,
        view = [55,250], alpha=1) 
plot_infls(valsf, ferrs, npatsf, fname='ferr', xlab = "$f_{err}$", line=False)
 
#floss, npatsfl, valsfl = scan_floss()
#pickle.dump([floss, npatsfl, valsfl], open('scan_floss_npats.pickled', 'wb'))  
floss, npatsfl, valsfl = pickle.load(open('scan_floss_npats.pickled', 'rb'))  
plot_3d(floss, npatsfl, valsfl, xlab="$f_{loss}$", ylab='patterns',
            zlab="$\epsilon$",fname='scan_floss_npat.png', vmin=0,
        view = [55,250], alpha=1)
plot_infls(valsfl, floss, npatsfl, fname='floss', xlab = "$f_{loss}$", line=False)

#Nheb, npatheb, valsheb = scan_vars(
#        Ns = [10, 300, 30], npats = [1, 100, 10],
#              fname = 'scan_N_npat_heb.png', method='heb')
    

#Nc, npatc, vals1, vals2 = test_learns()
#pickle.dump([Nc, npatc, vals1, vals2], open('scan_N_npats_hebopt.pickled', 'wb')) 
Nc, npatc, vals1, vals2 = pickle.load(open('scan_N_npats_hebopt.pickled', 'rb'))  
thresh = 0.138*Nc
plot_3d(Nc, npatc, vals1-vals2, xlab="N", ylab='patterns',
            zlab="$\epsilon_{heb} - \epsilon_{opt}$",fname='scan_N_npat_hebopt.png',
            vmin=0, alpha=1, zlim = [-0.05, 0.7])

plot_infls(vals1, Nc, npatc, vals2=vals2)
'''

#flosst, ferrst, valst = test_err_loss()
#pickle.dump([flosst, ferrst, valst], open('scan_floss_ferrs.pickled', 'wb')) 
flosst, ferrst, valst = pickle.load(open('scan_floss_ferrs.pickled', 'rb'))
plot_3d(flosst, ferrst, valst, xlab="$f_{loss}$", ylab='$f_{err}$',
            zlab="capacity",fname='test_err_loss.png', vmin=0,
            view = [45,70], alpha=1)
        
sts = ['hopfield', 'onomatopoeia', 'accommodate']
network = net(N = np.amax([len(st) for st in sts])*7)
network.learn_strings(sts)
#hop.query_string('neuroscence')
network.query_string('onomatopeia')
network.query_string('hobfield')
network.query_string('acommodate')

sts = ['i hope stephen gives me a good mark', 'i really like chocolate mousse',
       'i really like strawberries with milk']
network = net(N = np.amax([len(st) for st in sts])*7)
network.learn_strings(sts)
network.query_string('i hope stephen')
network.query_string('i really like')
network.query_string('strawberries with milk')
network.query_string('              strawberries with milk')




ferrso, npatsfo, valsfo = scan_ferrs(method = 'opt')
pickle.dump([ferrso, npatsfo, valsfo], open('scan_ferrs_npats_opt.pickled', 'wb'))
#ferrso, npatsfo, valsfo = pickle.load(open('scan_ferrs_npats_opt.pickled', 'rb'))  
plot_3d(ferrso, npatsfo, valsfo, xlab="$f_{err}$", ylab='patterns',
            zlab="$\epsilon$",fname='scan_ferr_npat_opt.png', vmin=0,
        view = [55,250], alpha=1) 
plot_infls(valsf, ferrs, npatsf, fname='ferr_opt', xlab = "$f_{err}$", line=False)


