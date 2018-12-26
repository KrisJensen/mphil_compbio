'''
After running screen_raw.py
we extract alignment counts and print and plot them
'''

import numpy as np
from parse_exonerate import main
import matplotlib.pyplot as plt

alignments = []
aligned_unaligned = []

for i in np.arange(70, 201, 1):
	print(i)
	files = ['BCc_genedict.pickled', 'exonerate_screen_raw/exonerate_BCc_coli_r'+str(i)+'.out']
	ind = [4,5,1,0]
	aligndict, res, alignees = main( [files, ind] )
	alignments.append(res['alignments'])
	aligned_unaligned.append(float(res['aligned'])**2 / float(res['alignments']) )
	print(aligned_unaligned[-1])

alignments = np.log10(np.array(alignments)+1)

plt.figure(figsize = (5,5))
plt.plot(np.arange(70, 201, 1), alignments )
plt.xlim([70, 200])
plt.ylim([0, max(alignments)])
plt.xlabel('Threshold (raw score)')
plt.ylabel('log10[alignments+1]')
plt.savefig('../figs/thresholding_raw_sq.png', dpi=560)
plt.show()

plt.figure(figsize = (6,4))
plt.plot(np.arange(70, 201, 1), aligned_unaligned )
plt.xlim([70, 200])
plt.ylim([0, 120])
plt.xlabel('Threshold (raw score)')
plt.ylabel('aligned^2 / alignments')
plt.savefig('../figs/thresholding_raw_frac_sq.png', dpi=560)
plt.show()



