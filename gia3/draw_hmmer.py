'''creates figure with visual display of aligned domains,
have boxes of different colours represent different domains
and have axis go from 1 to length of protein'''

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as patches

domains={'PF08934':'RB_c',\
	#strand-loop-helix: binds E2F1 and DP1
	'PF01858':'RB_A',\
	#cyclin fold of alpha helices
	'PF01857':'RB_B',\
	'PF11934':'DUF3452',\
	'PF00382':'TFIIB repeat'\
	}
cols = {'PF08934':'b',\
	'PF01858':'c',\
	'PF01857':'g',\
	'PF11934':'y',\
	'PF00382':'r'\
	}

#define proteins of interest and sort by similarity
specs = ['H. sapiens', 'C. jacchus', 'D. rerio', 'R. norvegicus', 'M. musculus',
	'G. gallus', 'F. catus', 'M. mulatta', 'I. punctatus', 'M. gallopavo',\
	'H. sapiens RBL1', 'H. sapiens RBL2']
sims = [100, 97.5, 53, 86.5, 90.6,
	72.4, 93.1, 96.8, 53.9, 71.2, 10, 5]
names = ['human', 'marmoset', 'zebrafish', 'rat', 'mouse',
	'chicken', 'cat', 'macaque', 'catfish', 'turkey', '', '']
specs = np.array(specs)
sims = np.array(sims)
args = np.argsort(-sims)
sims = sims[args]/100
specs = specs[args]
names = np.array(names)[args]


n = -1
k = 0
plt.figure(figsize = (6,4))

for line in open('hmmer/doms_parsed'):
	if line[0] == '>':#new protein
		n += 1
		if n%3 == 0:
			x = 0 #where to plot
			k -= 0.5 #change y coordinate
		elif n%3 == 1: x = 125 #second column
		elif n%3 == 2: x = 250 #third column

		l = int(line.split()[1][5:])
		if n == 0: l0 = l #define H.sapiens RB1 to have length 100
		plt.plot( [0+x, 100*l/l0+x], [k, k] , 'k-', linewidth=2) #plot length of protein
		plt.text(x+50, k-0.17, specs[n],\
			horizontalalignment = 'center',\
			verticalalignment = 'center',\
			fontstyle = 'italic')
		plt.text(x+50, k-0.27, names[n],\
			horizontalalignment = 'center',\
			verticalalignment = 'center')#,\
			#fontstyle = 'italic')

	else: #plot a box corresponding to this domain
		start, end = [int(p) for p in line.split()[1].split(':')[1].split(',')]
		plt.fill( [100*start/l0+x, 100*end/l0+x, 100*end/l0+x, 100*start/l0+x],
				[k-0.1, k-0.1, k+0.1, k+0.1], cols[line.split()[0].split(':')[1]]+'-' )
#save figure
plt.axis('off')
plt.savefig('/local/data/public/ktj21/GIA3/figures/domains.png', dpi = 360, bbox_inches = 'tight')
plt.show()
plt.close()
