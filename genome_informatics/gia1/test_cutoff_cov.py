'''
After settling on a hash length, we want to scan the cutoff/exp_cov landscape. We do this by trying a number off different combinations of values
'''

import numpy as np
from subprocess import call

def try_cut_cov(directory, cutoffs = np.arange(0, 25.5, 0.5), covs = np.arange(0, 51, 1)):
	'''given the name of a hashed directory, apply a number of different coverages and cutoffs'''

	for cut in cutoffs: #hardcoded; from no cutoff to our max coverage
		for cov in covs: #exp_cov. Also goes from zero to beyond our max coverage
			call(['/local/data/public/ktj21/programs/velvet/velvetg', directory, '-cov_cutoff', str(cut), '-exp_cov', str(cov), '-ins_length', '500']) #run velvetg

if __name__ == '__main__':
	try_cut_cov('cut_cov_57')

