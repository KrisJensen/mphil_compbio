"""
code for generating figure 1 in Cancer Evolution assignment 1.
Additional annotations were added in LaTeX.
The graph (figure 1a) was visualized in cytoscape.
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
import networkx as nx
import copy

G = nx.Graph() #initialize graph

#prepare figure 1b
ax = plt.subplots()[1]
plt.xlim(0,40)
plt.ylim(0, 25)

#initialize vectors for storing data for summary plots
fracs = []
xs = []
nmuts = []
muts = []
names = []
ncells= []
cols = ['#d3d3d3', 'g', 'b', 'y', 'c', 'm', '#00FF00', '#9933FF', '#FF66FF' ]

def rects(lims, col, name, parent, size):
    '''given x, y values for a clonal block, plots it and stores statistics'''
    rect = Rectangle((lims[0], lims[2]), lims[1]-lims[0],\
                     lims[3]-lims[2], facecolor = col)
    ax.add_patch(rect)
            
    #store some information for row/column summaries
    fracs.append((lims[3]-lims[2])/25)
    xs.append([lims[0], lims[1]])
    nmuts.append(lims[1]-lims[0])
    names.append(name)
    muts.append(name[-1])
    ncells.append(size)
    
    #add node to graph
    G.add_node(name, lab = str(int(size/25*100)), size = int(size/25*100)) 
    if not parent == 'None':
        #add edge to graph
        G.add_edge(parent, name, label = '+'+name[-1])
            
#add segments corresponding to each clone
rects([0,40,0,25], '#d3d3d3', '0', 'None', 5)  
rects([0,5,0,20], 'g', 'A', '0', 5)    
rects([5,5+11,0,15], 'b', 'AB', 'A', 5)
rects([16,16+3,0,10], 'y', 'ABC', 'AB', 0)     
rects([19,19+6,0,6], 'c', 'ABCD', 'ABC', 0)
rects([25,25+6,6,6+4], 'm', 'ABCE', 'ABC', 1)
rects([31,31+1,6,6+3], '#00FF00', 'ABCEF', 'ABCE', 3)
rects([32,32+5,0,3], '#9933FF', 'ABCDG', 'ABCD', 3)   
rects([37,37+3,3,3+3], '#FF66FF', 'ABCDH', 'ABCD', 3)
      
#add gridlines
for x in range(40):
    plt.axvline(x, color='w')
for y in range(25):
    plt.axhline(y, color='w')

#save figure
dirn = '/Users/kris/Documents/cambridge/mphil/assignments/ce/ce1/'
plt.xticks([])
plt.yticks([])
plt.savefig(dirn+'squares.png', bbox_inches = 'tight', pad_inches = 0.8)
plt.show()

#store graph
nx.write_graphml(G,
    '/Users/kris/Documents/cambridge/mphil/assignments/ce/ce1/graph.graphml')
nx.draw(G, with_labels = True)
plt.show()

fs = 16 #fontsize in plots
s = 0.8 #scaling of plot size

#create row summary
inds = np.argsort(ncells)
plt.figure(figsize = np.array([6,3])*s)
plt.bar(range(len(names)), np.array(ncells)[inds],\
        color = np.array(cols)[inds], width = 0.2)
plt.xticks(range(len(names)), np.array(names)[inds], rotation=90, weight='bold')
plt.yticks(FontSize = fs-4)
plt.title('Row summary', FontSize=fs)
plt.xlabel('Cellular populations', FontSize=fs)
plt.ylabel('Number of cells (25$\pi$)', FontSize=fs)
plt.savefig(dirn+'rows.png', bbox_inches = 'tight', dpi = 360)
plt.show()

#create column summary
inds = np.argsort(-np.array(nmuts[1:]))
newfracs = copy.deepcopy(fracs); newfracs[5]=0.18; newfracs[4] = 0.255
ticks, va = [], []
for f in np.array(fracs[1:])[inds]: #need ticks for cellular fraction and CCF
    ticks.append(str(int(100*f)))
    va.append(-0.03)
    ticks.append('('+str(int(100*f*25/20))+')')
    va.append(-0.14)
    ticks.append(' ')
    va.append(-0.23)
fig = plt.figure(figsize = np.array([6,3])*s )
ax = plt.gca()
plt.bar(np.array([int(f*100) for f in newfracs[1:]])[inds],
        np.array(nmuts[1:])[inds],
        color = np.array(cols[1:])[inds], width = 1.5)
ax.set_xticks( np.repeat(np.array([int(f*100) for f in newfracs[1:]])[inds], 3))
ax.set_xticklabels( ticks )
for t, y in zip( ax.get_xticklabels( ), va ):
    t.set_y( y )
plt.yticks(FontSize = fs-4)
plt.title('Column summary', FontSize=fs)
plt.xlabel('Cellular fraction (CCF)', FontSize=fs)
plt.ylabel('# mutations', FontSize=fs)
plt.savefig(dirn+'cols.png', bbox_inches = 'tight', dpi=360)
plt.show()

#create 'bulk sequencing' graph
nmuts = [9,6,6,3,11,5]
fracs= [0.12,0.16,0.24,0.40,0.60,0.80]
ticks, va = [], []
for f in np.array(fracs):
    ticks.append(str(int(100*f)))
    va.append(-0.02)
    ticks.append(' ')
    va.append(-0.02)
    ticks.append(' ')
    va.append(-0.0)
fig = plt.figure(figsize = np.array([4,3])*s*1.3 )
ax = plt.gca()
plt.bar(np.array([int(f*100) for f in fracs]),
        np.array(nmuts),
        color = '#d3d3d3', width = 2)
ax.set_xticks( np.repeat(np.array([int(f*100) for f in fracs]), 3))
ax.set_xticklabels( ticks )
for t, y in zip( ax.get_xticklabels( ), va ):
    t.set_y( y )
plt.yticks(FontSize = fs-4)
plt.title('Column summary', FontSize=fs)
plt.xlabel('Cellular fraction', FontSize=fs)
plt.ylabel('# mutations', FontSize=fs)
plt.savefig(dirn+'bulk.png', bbox_inches = 'tight', dpi=360)
plt.show()


