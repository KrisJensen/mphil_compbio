'''
Read logs from screening of kmer lengths and cutoffs to plot N50 as a function of both parameters and return the optimum value
'''

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm

if __name__ == '__main__':
	N50_all = [] #matrix of N50 values
	cutoffs_all = []
	kmer_all = []
	for kmer in range(17, 97, 2): #same as other scripts
		f = open('kmer_'+str(kmer)+'/Log', 'r') #read log
		cutoffs, N50s = [], []
		for line in f:
			split =line.split()
			if len(split) > 0 and 'velvetg' in split[0]: #this line specifies parameters including cutoff
				cutoffs.append( float(split[3]) )
			if len(split) > 0 and split[0] == 'Final': #this line specifies result including N50
				N50s.append( int(split[8][0:-1]) )

                #N50_all += N50s
		#cutoffs_all += cutoffs
		#kmer_all += [kmer for x in range(len(N50s))]

		N50_all.append(N50s) #add row of values for given kmer length to matrix
		cutoffs_all.append(cutoffs)
		kmer_all.append([kmer for x in range(len(N50s))]) #create grid for plt plotting
		f.close() #ready for next kmer length

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d') #3d plot

	kmer_all = np.array(kmer_all)
	N50_all = np.log10(np.array(N50_all)) #plot log10(N50) as we have order of magnitude variation
	cutoffs_all = np.array(cutoffs_all)

	print(N50_all.shape, kmer_all.shape, cutoffs_all.shape)

	maxN = np.amax(N50_all) #find our optimum value
	print('max N50 is', np.round(10**maxN, 0), '\ncutoff', cutoffs_all[ N50_all == maxN ], '\nkmer', kmer_all[ N50_all == maxN ])

        #ax.scatter(kmer_all, cutoffs_all, np.log10(N50_all))
	surf = ax.plot_surface(kmer_all, cutoffs_all, N50_all, cmap=cm.coolwarm, antialiased = False, vmin = 0, vmax = np.amax(N50_all))
	fig.colorbar(surf, shrink=0.5, aspect=6)

	#make things look nice
	ax.set_xlabel('kmer length / bp')
	ax.set_ylabel('cov cutoff')
	ax.set_zlabel('log10(N50)')
	plt.xlim(np.amin(kmer_all), np.amax(kmer_all))
	plt.ylim(np.amin(cutoffs_all), np.amax(cutoffs_all))

	#save and show
	plt.savefig('test_surface_kmer_cutoff.png')
	plt.show()          



