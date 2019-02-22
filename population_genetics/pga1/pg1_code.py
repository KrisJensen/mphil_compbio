#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 20:00:55 2019

@author: kris
"""

import matplotlib
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from scipy.misc import factorial


t = np.array(range(7))

n1 = np.array([28, 39, 57, 72, 78, 81, 82])
n0 = 100 - n1

N = 100

sigs = np.linspace(0, 2, N)
q0s = np.linspace(0.05,0.93,N)


prefactor = factorial(100) / (factorial(n1) * factorial(n0))

Ls = np.zeros((N,N))

def calcl(sig, q0, gamma = 1, t=t, prefactor=prefactor, n1=n1, n0=n0, ms = False):
    #q = (q0*np.exp(sig * t)) / (1 - q0 + q0*np.exp(sig * t))
    q = gamma * (q0*np.exp(sig * gamma * t)) / (gamma - q0 + q0*np.exp(sig *gamma* t))
    L = np.sum( np.log( prefactor * q**n1 * (1-q)**n0 ) )
    if ms:
        msq = (np.mean((100*q-n1)**2))**0.5
        print(msq)
        print(np.log( prefactor * q**n1 * (1-q)**n0 ))
    return L

for i, sig in enumerate(sigs):
    for j, q0 in enumerate(q0s):
        Ls[i,j] = calcl(sig, q0)

sigs = np.repeat(sigs,N).reshape(N,N)
q0s = np.transpose(np.repeat(q0s,N).reshape(N,N))

fig = plt.figure(figsize = (7,4))
ax = fig.add_subplot(111, projection='3d')

Ls = -Ls
surf = ax.plot_surface(sigs, q0s, Ls, cmap=cm.Greys, antialiased = False,
                       vmin = 0, vmax = np.amax(Ls))
fig.colorbar(surf, shrink=0.5, aspect=6, pad=0.05)

#make things look nice
ax.set_xlabel('$\sigma$')
ax.set_ylabel('$q_0$')
ax.set_zlabel('$-L$')
plt.xlim(0,2)
plt.ylim(0,1)
#save and show
ax.view_init(45, 325)
ax.dist = 10.5
plt.savefig('Lsurf.png', dpi=360)

plt.show()


def sdescent(q0 = 0.05, sig = 0, gamma=1, d = 10**(-5), g = 10**(-7),
             thresh = 10**(-12), pas = False):
    
    
    L = -calcl(sig, q0, gamma=gamma)
    gammas = [gamma]
    
    err = thresh+1
    
    sigs = [sig]
    q0s = [q0]
    Ls = [L]
    
    niter = 0
    print(niter, sig, q0, L)
    while err > thresh and niter < 10**9:
        niter += 1
        
        L0 = L
        
        newsig = sig+d
        newq = q0+d
        

        d1 = -calcl(newsig, q0, gamma=gamma)
        d2 = -calcl(sig, newq, gamma=gamma)
            
        if pas:
            newgam = gamma+d
            d3 = -calcl(sig, q0, gamma=newgam)
            gamma = gamma + (L - d3)/d*g
            gammas.append(gamma)
        
        sig = sig + (L - d1)/d*g
        q0 = q0 + (L - d2)/d*g

        L = -calcl(sig, q0, gamma=gamma)
        
        sigs.append(sig)
        q0s.append(q0)
        Ls.append(L)
        
        err = np.abs(L0-L)
        if niter % 10000 == 0: print(niter, sig, q0, gamma, L, err)
            

    if pas:    
        print(sig, q0, gamma, L)    
        return sig, q0, gamma, L, sigs, q0s, gammas, Ls
    
    else:
        print(sig, q0, L)    
        return sig, q0, L, sigs, q0s, Ls
    
#sigmin, q0min, Lmin, sigsp, q0sp, Lsp = sdescent()
#sigmint, q0mint, gmint, Lmint, sigspt, q0spt, gspt, Lspt = sdescent(pas = True)

def fit_model(sig = sigmin, q0=q0min, gamma=1.0, ts = np.linspace(0,7,100),
              ext = '', pas=False):
    
    if pas:n = 100 * gamma * \
        (q0*np.exp(sig * gamma * ts)) / (gamma - q0 + q0*np.exp(sig *gamma* ts))
    else: n = (q0*np.exp(sig * ts)) / (1 - q0 + q0*np.exp(sig * ts))*100
    plt.plot(t, n1/100, 'bx')
    plt.plot(ts, n/100, 'y--')
    plt.xlabel('time/ days')
    plt.ylabel('$q_i^1$')
    plt.legend(['Data','Model'])
    plt.savefig('fit_model'+ext+'.png', dpi=360)
    plt.show()

fit_model()

plt.plot(t, n1, 'bx')
plt.xlabel('time/ days')
plt.ylabel('$n_i^1$')
plt.savefig('data.png', dpi=360)
plt.show()

######################################

fig = plt.figure(figsize = (7,4))
ax = fig.add_subplot(111, projection='3d')

surf = ax.plot_surface(sigs, q0s, Ls, cmap=cm.Greys, antialiased = False,
                       vmin = 0, vmax = np.amax(Ls))
fig.colorbar(surf, shrink=0.5, aspect=6, pad=0.05)

#make things look nice
ax.set_xlabel('$\sigma$')
ax.set_ylabel('$q_0$')
ax.set_zlabel('$-L$')
plt.xlim(0,2)
plt.ylim(0,1)
#save and show
ax.view_init(45, 355)
ax.dist = 10.5

plt.plot(sigsp, q0sp, np.array(Lsp)+1, 'k-')

plt.savefig('Lsurfproj.png', dpi=360)

plt.show()

######################################

fig = plt.figure(figsize = (7,4))
ax = fig.add_subplot(111, projection='3d')

surf = ax.plot_surface(sigs, q0s, Ls, cmap=cm.Greys, antialiased = False,
                       vmin = 0, vmax = np.amax(Ls))
fig.colorbar(surf, shrink=0.5, aspect=6, pad=0.05)

#make things look nice
ax.set_xlabel('$\sigma$')
ax.set_ylabel('$q_0$')
ax.set_zlabel('$-L$')
plt.xlim(0,2)
plt.ylim(0,1)
#save and show
ax.view_init(45, 355)
ax.dist = 10.5

Lthresh = Lmin+np.log(100)
coords = (Ls < Lthresh) #& (Ls > Lthresh-0.5)
sthresh = sigs[coords]
q0thresh = q0s[coords]
Lsthresh = Ls[coords]

print(np.amin(sthresh), sigmin, np.amax(sthresh))
print(np.amin(q0thresh), q0min, np.amax(q0thresh))

plt.plot(sthresh, q0thresh, Lsthresh, 'k-')

plt.savefig('Lsurfthresh.png', dpi=360)

plt.show()

#########

l0 = -calcl(sigmin, q0min)
delta = 0.001
ls = -calcl(sigmin+delta, q0min)
lq = -calcl(sigmin, q0min+delta)

print((ls-l0)/delta, (lq-l0)/delta)

svec = np.array([(ls-l0)/delta, (lq-l0)/delta])
tvec = np.array([(l0-lq)/delta, (ls-l0)/delta])
svec = svec/np.linalg.norm(svec); tvec = tvec/np.linalg.norm(tvec)
mat = np.array([svec, tvec])

np.dot(mat, np.array([0,0.05]))
np.dot(mat, np.array([2,0.05]))
np.dot(mat, np.array([0,0.93]))
np.dot(mat, np.array([2,0.93]))

minv = np.linalg.inv(mat)


N = 100

scoords = np.linspace(0.29,0.85, N)
tcoords = np.linspace(-0.552,0.15,N)


prefactor = factorial(100) / (factorial(n1) * factorial(n0))

Lnews = np.zeros((N,N))

for i, scoord in enumerate(scoords):
    for j, tcoord in enumerate(tcoords):
        sig, q0 = np.dot(minv, np.array([scoord,tcoord])) 
        #print(scoord, tcoord, sig, q0)
        Lnews[i,j] = calcl(sig, q0)

ss_p = np.repeat(scoords,N).reshape(N,N)
ts_p = np.transpose(np.repeat(tcoords,N).reshape(N,N))

fig = plt.figure(figsize = (7,4))
ax = fig.add_subplot(111, projection='3d')

Lnews = -Lnews
surf = ax.plot_surface(ss_p, ts_p, Lnews, cmap=cm.Greys, antialiased = False,
                       vmin = 0, vmax = np.amax(Lnews))
fig.colorbar(surf, shrink=0.5, aspect=6, pad=0.05)

smin, tmin = np.dot(mat, np.array([sigmin, q0min]))

plt.plot([smin], [tmin], [Lmin], 'ko')

coords = (Lnews < Lthresh) #& (Ls > Lthresh-0.5)
ssthresh = ss_p[coords]
ttthresh = ts_p[coords]
LLsthresh = Lnews[coords]

print(np.amin(ssthresh), smin, np.amax(ssthresh))
print(np.amin(ttthresh), tmin, np.amax(ttthresh))

plt.plot(ssthresh, ttthresh, LLsthresh, 'k-')

#make things look nice
ax.set_xlabel('$s$')
ax.set_ylabel('$t$')
ax.set_zlabel('$-L$')
plt.xlim(0.29,0.85)
plt.ylim(-0.552,0.15)
#save and show
ax.view_init(45, 45)
ax.dist = 10.5
plt.savefig('Lsurfrot.png', dpi=360)

plt.show()

sig1, q01 = np.dot(minv, np.array([smin, tmin-0.16]))
fit_model(sig=sig1, q0=q01, ext = '1')
sig2, q02 = np.dot(minv, np.array([smin, tmin+0.14]))
fit_model(sig=sig2, q0=q02, ext = '2')
sig3, q03 = np.dot(minv, np.array([smin, tmin-0.32]))
fit_model(sig=sig3, q0=q03, ext = '3')

fit_model(sig=sigmint, q0=q0mint, gamma=gmint, pas=True, ext = 'pas')


#############################

sigst = np.linspace(0, 2, N)
q0st = np.linspace(0.12,0.93,N)

Lst = np.zeros((N,N))

for i, sig in enumerate(sigst):
    for j, q0 in enumerate(q0st):
        Lst[i,j] = calcl(sig, q0, gamma = gmint)

sigst = np.repeat(sigst,N).reshape(N,N)
q0st = np.transpose(np.repeat(q0st,N).reshape(N,N))

fig = plt.figure(figsize = (7,4))
ax = fig.add_subplot(111, projection='3d')

Lst = -Lst


Lthresht = Lmint+np.log(100)
coordst = (Lst < Lthresht) #& (Ls > Lthresh-0.5)
sthresht = sigst[coordst]
q0thresht = q0st[coordst]
Lsthresht = Lst[coordst]

print(np.amin(sthresht), sigmint, np.amax(sthresht))
print(np.amin(q0thresht), q0mint, np.amax(q0thresht))

surf = ax.plot_surface(sigst, q0st, Lst, cmap=cm.Greys, antialiased = False,
                       vmin = 0, vmax = np.amax(Lst))
fig.colorbar(surf, shrink=0.5, aspect=6, pad=0.05)
plt.plot(sthresht, q0thresht, Lsthresht, 'k-')
#make things look nice
ax.set_xlabel('$\sigma$')
ax.set_ylabel('$q_0$')
ax.set_zlabel('$-L$')
plt.xlim(0,2)
plt.ylim(0.1,1)
#save and show
ax.view_init(45, 25)
ax.dist = 10.5
plt.savefig('Lsurfpas.png', dpi=360)

plt.show()


##################

'''
sigst = np.linspace(0, 2, N)
q0st = np.linspace(0.12,0.93,N)
gst = np.linspace(0.15,1.0,N)

Lst = np.zeros((N,N,N))

for i, sig in enumerate(sigst):
    print('new i', i)
    for j, q0 in enumerate(q0st):
        for k, gamma in enumerate(gst):
            Lst[i,j,k] = calcl(sig, q0, gamma = gamma)
            
print('filled out')

sigst =  np.repeat(sigst,N*N).reshape(N,N,N)
q0st = np.repeat( np.transpose(np.repeat(q0st,N).reshape(N,N)), N).reshape(N,N,N)
gst =  np.transpose(np.repeat(gst,N**2).reshape(N,N,N))

print('made arrays')

Lst = -Lst


Lthresht = Lmint+np.log(100)
coordst = (Lst < Lthresht) #& (Ls > Lthresh-0.5)
sthresht = sigst[coordst]
q0thresht = q0st[coordst]
gthresht = gst[coordst]
Lsthresht = Lst[coordst]

print('found threshs')

print(np.amin(sthresht), sigmint, np.amax(sthresht))
print(np.amin(q0thresht), q0mint, np.amax(q0thresht))
print(np.amin(gthresht), gmint, np.amax(gthresht))
'''


#######################################


dirs = ['1_mut', '2_mut', '3_mut']
newlines = [ [87], [1, 97], [1,74,112], [], []  ]
cols = [ ['b'], ['y', 'b'], ['y', 'b', 'b'], [], [] ]

for n in range(1,6):

    f1 = open('/Applications/optimist/'+str(n)+'_mut/Model_frequencies.out', 'r')
    f2 = open('/Applications/optimist/'+str(n)+'_mut/Real_frequencies.out', 'r')
    
    dat = []
    ltypes = ['-', 'x']
    for i, f in enumerate([f1, f2]):
        ts = []
        fr1 = []
        fr2 = []
        for line in f:
            split = [float(val) for val in line.split()]
            ts.append(split[0])
            fr1.append(split[1])
            fr2.append(split[2])
    
        dat.append([ts, fr1, fr2])
    
        plt.plot(ts, fr1, 'b'+ltypes[i] )
        plt.plot(ts, fr2, 'y'+ltypes[i] )
    
    for i, val in enumerate(newlines[n-1]):
        plt.axvline(x=val, color = cols[n-1][i], linestyle = ':')
        
    plt.xlabel('generation')
    plt.ylabel('frequency')
    plt.legend(['Model frequency 1', 'Model frequency 2', 'Real frequency 1',
                'Real frequency2'])
    plt.savefig('fit_mut'+str(n)+'.png', dpi=360)
    plt.show()
    
    f1.close()
    f2.close()
    
muts = []
lls = []
with open('/Applications/optimist/systematic_lls', 'r') as f:
    for i, line in enumerate(f):
        muts.append(i+1)
        lls.append(float(line))

AIC = 2.5 * 2*np.array(muts) - np.array(lls)
plt.plot(muts, -np.array(lls), 'b-')
plt.plot(muts, AIC, 'y-')
plt.legend(['-(Log likelihood)', 'AIC'])
plt.xlabel('number of mutations')
plt.ylabel('likelihood score')
plt.savefig('lls.png', dpi = 360)
plt.show()



pi = 0.9
a = 0.8
b = 0.9
c = 0.7

p = np.linspace(0,1,1000)

f = pi*(p*p*a + 2*(1-p)*p*b + (1-p)*(1-p)*c) + (1-pi)*(p*p + 2*(1-p)*p + (1-p)*(1-p)*c)

plt.plot(p, f, 'b-')
plt.xlim([0,1])
plt.xlabel('p')
plt.ylabel('mean fitness')
#plt.savefig('fitness.png', dpi=360)
plt.show()

df=pi*(2*p*a + 2*(1-2*p)*b + (2*p-2)*c)+(1-pi)* (2*p + 2*(1-2*p) + (2*p-2)*c)

plt.plot(p, df, 'b-')
plt.plot(p,np.repeat(0.0, len(p)), 'k-')
plt.xlim([0,1])
plt.xlabel('p')
plt.ylabel('dF')
#plt.savefig('dfitness.png', dpi=360)
plt.show()

N = 1000

b = np.linspace(0,1,N)
b = np.repeat(b,N).reshape(N,N)
p = np.linspace(0,1,1000)
p = np.transpose(np.repeat(p,N).reshape(N,N))

f = pi*(p*p*a + 2*(1-p)*p*b + (1-p)*(1-p)*c) + (1-pi)*(p*p + 2*(1-p)*p + (1-p)*(1-p)*c)

fig = plt.figure(figsize = (7,4))
ax = fig.add_subplot(111, projection='3d')

surf = ax.plot_surface(b, p, f, cmap=cm.Greys, antialiased = False,
                       vmin = np.amin(f), vmax = np.amax(f))
fig.colorbar(surf, shrink=0.5, aspect=6, pad=0.05)

#make things look nice
ax.set_xlabel('b')
ax.set_ylabel('p')
ax.set_zlabel('fitness')
plt.xlim(0,1)
plt.ylim(0,1)
#save and show
ax.view_init(45, 325)
ax.dist = 10.5
plt.savefig('fsurf.png', dpi=720)

plt.show()

cols = np.zeros((N, N))

for i in range(N):
    for j in range(N):
        bc = b[i,0]
        pc = p[0,j]
        
        if bc < 2/3:
            if pc > (1.2 - 1.8*bc)/(2.64 - 3.6*bc): cols[i,j] = 1.0
            else: cols[i,j] = 0.0
            
        elif bc > 0.8:
            cols[i,j] = (1.2 - 1.8*bc)/(2.64 - 3.6*bc)
            
        else:
            cols[i,j] = 1.0
            
print('finished cols')

newfig = plt.figure()
newax = newfig.add_subplot(111, projection='3d')
surf1 = newax.plot_surface(b, p, cols, cmap=cm.Greys, antialiased = False,
                           vmin = 0, vmax = 1)
plt.show()       

minn, maxx = cols.min(), cols.max()
norm = matplotlib.colors.Normalize(minn, maxx)
m = plt.cm.ScalarMappable(norm=norm, cmap=cm.Greys)
m.set_array([])
fcolors = m.to_rgba(cols)

fig = plt.figure(figsize = (7,4))
ax = fig.add_subplot(111, projection='3d')
print('about to plot')
surf = ax.plot_surface(b,p,f,facecolors=fcolors, vmin=minn, vmax=maxx, shade=False)
fig.colorbar(surf1, shrink=0.5, aspect=6, pad=0.05)
ax.view_init(45, 325)
ax.dist = 10.5


#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#surf = ax.plot_surface(b, p, f, cmap=cm.Greys, antialiased = False,
#    vmin = np.amin(f), vmax = np.amax(f))

#make things look nice
ax.set_xlabel('b')
ax.set_ylabel('p')
#ax.set_zlabel('final composition')
plt.xlim(0,1)
plt.ylim(0,1)
#save and show
print('saving')
plt.savefig('fsurfinit.png', dpi=720, bbox_inches='tight')
plt.show()






