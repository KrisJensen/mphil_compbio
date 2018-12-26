
'''considers variation of the protein from ensembl data.avg of 4 mutations per residue; do bins of 3 amino acids'''

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

binsize = 3
rcounts = [0 for i in range(929)]
counts = [0 for i in range(int(929/binsize)+1)]
with open('RB1.pep.var.csv', 'r') as f:

	for line in f:
		try:
			split = line.split(',')
			if split[3] == '"coding sequence variant"':	
				counts[ int(float(split[0])/binsize) ] += 1
				rcounts[ int(split[0])-1 ] += 1
		except:
			print('error\n', line)
counts = np.array(counts)
rcounts = np.array(rcounts)

#phos = np.array([ 249, 252, 356, 373, 608, 780, 788, 795, 807, 826]) #these are the ones that are fully conserved
#phos = np.array([230, 249, 252, 356, 373, 608, 612, 780, 788, 795, 807, 811, 821, 826])
phos = np.array([249, 252, 356, 373, 608, 612, 780, 788, 795, 807, 811, 821, 826])

bphos = np.array([int(i/binsize) for i in phos])
phos = phos-1
print(rcounts[phos])
pcounts = rcounts[phos]

plt.plot(counts, 'bo', markersize = 1.5)
plt.plot(bphos, counts[bphos], 'ro', markersize = 2)
plt.show()
plt.close()

plt.figure(figsize = (5,3.5))
plt.plot(rcounts, 'bo', markersize = 1.5)
plt.plot(phos, pcounts, 'ro', markersize = 2)
plt.xlabel('residue')
plt.ylabel('ensembl coding sequence mutations')
plt.savefig('/local/data/public/ktj21/GIA3/figures/pepvar_distribution.png', dpi=360, bbox_inches = 'tight')
plt.show()

plt.close()

plt.figure(figsize = (5,3.5))
bins = range(max(rcounts)+1)
plt.hist(rcounts, bins = bins, density=True, color=[0,0,1,1] )
plt.hist(pcounts, bins = bins, density=True, color = [0,1,0,0.5] )
plt.legend(['Other residues', 'Phosphorylation sites'])
plt.xlabel('ensembl coding sequence mutations')
plt.ylabel('frequency')
plt.savefig('/local/data/public/ktj21/GIA3/figures/pepvar_hist.png', dpi=360, bbox_inches='tight')
plt.show()
print(np.mean(rcounts), np.mean(rcounts[phos] ))
print(stats.ttest_ind(rcounts, pcounts, equal_var = False) ) #welch's t-test
print(' '.join(['' for i in range(16)])+'  '.join([ str(c) for c in pcounts]))



