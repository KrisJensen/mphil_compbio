#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
simple consideration of number of ancestors
"""

import numpy as np
import matplotlib.pyplot as plt

gs = range(1,16)
ns = [17*(2**g) for g in gs]

ps = []
N = 1.2*10**8
N /= 3 #generation is third of population

for n in ns:
    P = 1
    for i in range(n):
        P = P * (N-i)/N
    ps.append(1-P)
    
plt.figure(figsize = (5,3))
plt.plot(gs, ps, 'b-')
plt.xlabel('time / generations')
plt.ylabel('$P_{shared}$')
plt.savefig('figures/simple_ancestors.png', dpi=120, bbox_inches='tight')
plt.show()

for i in range(len(gs)):
    print(gs[i], gs[i]*25, ps[i])
