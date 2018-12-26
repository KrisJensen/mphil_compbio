'''
Script for reading logs of screen of exp_cov and cutoff values and plotting N50 as a function of these parameters.
Similar to read_logs.py, but the way we parse the data is different enough to warrant two scripts rather than a long list of conditional statements in my opinion.
'''

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm

if __name__ == '__main__':
	N50_all = [] #matrix of N50 values
	cutoffs_all = []
	covs_all = []
	f = open('kmer31_cut_cov_screen/Log', 'r') #read log
	covs, N50s, cutoffs = [], [], []

	cov0, cut0 = -1, -1

	for line in f:
		split =line.split()
		if len(split) > 0 and split[0] == 'velvetg': #specify parameters in this line
			cutoff, cov = float(split[3]), float(split[5])
		if len(split) > 0 and split[0] == 'Final': #get result in this line
			N50 = int(split[8][0:-1])

			if cov < cov0: #if we go back to cov=0, we have reached a new row of the matrix and need to process the previous row
				cov0 = -1 #reset counter
				N50_all.append(N50s) #add row to matrix
				cutoffs_all.append(cutoffs)
				covs_all.append(covs)
				N50s, cutoffs, covs = [N50], [cutoff], [cov] #initialise new row
			else: #still in same row so just add values
				cutoffs.append(cutoff)
				covs.append(cov)
				N50s.append(N50)
				cov0, cut0 = cov, cutoff
	f.close()

               #N50_all += N50s
	#cutoffs_all += cutoffs
	#kmer_all += [kmer for x in range(len(N50s))]

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')

	covs_all = np.array(covs_all)
	N50_all = np.log10(np.array(N50_all)) #plot log10 as we ahve order or magnitude variation
	cutoffs_all = np.array(cutoffs_all)
	
	maxN = np.amax(N50_all) #find our optimum value
	print('max N50 is', np.round(10**maxN, 0), '\ncutoff', cutoffs_all[ N50_all == maxN ], '\nexp_cov', covs_all[ N50_all == maxN ])          

	#ax.scatter(kmer_all, cutoffs_all, np.log10(N50_all))
	surf = ax.plot_surface(covs_all, cutoffs_all, N50_all, cmap=cm.coolwarm, antialiased = False,vmin = 0, vmax = np.amax(N50_all))
	fig.colorbar(surf, shrink=0.5, aspect=6)

	#make things look nice
	ax.set_xlabel('exp_cov')
	ax.set_ylabel('cov_cutoff')
	ax.set_zlabel('log10(N50)')
	plt.xlim(np.amin(covs_all), np.amax(covs_all))
	plt.ylim(np.amin(cutoffs_all), np.amax(cutoffs_all))

	#save
	plt.savefig('test_surface_cov_cutoff.png')
	plt.show()




