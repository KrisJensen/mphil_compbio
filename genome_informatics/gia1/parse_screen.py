'''
After running screen_exonerate.py
we extract alignment counts and print and plot them
'''

import numpy as np
from parse_exonerate import main
import matplotlib.pyplot as plt

alignments = []
aligned_unaligned = []

for i in np.arange(0.5, 100.5, 0.5):
	print(i)
	files = ['BCc_genedict.pickled', 'exonerate_screen_2/exonerate_BCc_coli_p'+str(i)+'.out']
	ind = [4,5,1,0]
	aligndict, res, alignees = main( [files, ind] )
	alignments.append(res['alignments'])
	if res['alignments'] == 0: aligned_unaligned.append(0)
	else: aligned_unaligned.append( float(res['aligned'])**2 / float(res['alignments']) )
	print(aligned_unaligned[-1])

alignments = np.log10(np.array(alignments)+1)

plt.figure(figsize = (6,4))
plt.plot(np.arange(0.5, 100.5, 0.5), alignments )
plt.xlim([0, 60])
plt.ylim([0, max(alignments)])
plt.xlabel('Threshold (percent max score)')
plt.ylabel('log10[alignments+1]')
plt.savefig('../figs/thresholding_sq.png', dpi=540)
plt.show()



plt.figure(figsize = (5,5))
plt.plot(np.arange(0.5, 100.5, 0.5), aligned_unaligned )
plt.xlim([0, 60])
plt.ylim([0, 100])
plt.xlabel('Threshold (raw score)')
plt.ylabel('aligned^2 / alignments')
plt.savefig('../figs/thresholding_p_frac_sq.png', dpi=540)
plt.show()



