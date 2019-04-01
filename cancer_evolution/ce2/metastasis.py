#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Code for generating CCF plots for metastatic trees
"""

import matplotlib.pyplot as plt

#CCFs for reference tree
ccfs01 = [100, 75, 50, 30, 20, 15, 15, 15, 0, 0, 0 ]
ccfs02 = [100, 75, 50, 30, 20, 15, 15, 15, 0, 0 ]

#clonal colours
cols1 = ['g', 'b', 'y', 'c', 'm', '#00FF00', '#9933FF', '#FF66FF', 'r', 'r', 'r' ]
cols2 = ['g', 'b', 'y', 'c', 'm', '#00FF00', '#9933FF', '#FF66FF', 'r', 'r' ]
    
#CCFs for metastatic tumors
ccfs1 = [100, 0 , 0 , 0 , 0 , 0 , 0 , 0 , 75, 30, 20]
ccfs2 = [100, 100, 100, 0, 100, 0, 0, 0, 60, 30]

#tree 1
plt.figure(figsize = (3,3))
plt.scatter(ccfs01, ccfs1, color=cols1, alpha=0.8, s = 50)
plt.xlabel('primary')
plt.ylabel('metastasis')
plt.savefig('meta_early.png', dpi=120, bbox_inches = 'tight')
plt.show()

#tree 2
plt.figure(figsize = (3,3))
plt.scatter(ccfs02, ccfs2, color=cols2, alpha=0.8, s=50)
plt.xlabel('primary')
plt.ylabel('metastasis')
plt.savefig('meta_late.png', dpi=120, bbox_inches = 'tight')
plt.show()

