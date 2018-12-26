import matplotlib.pyplot as plt
import numpy as np
from scipy.misc import comb
import sys

BL50 = {} #construct BLOSUM50 matrix
with open('/local/data/public/ktj21/programs/fasta36/data/blosum50.mat', 'r') as f:
	for i in range(5): f.readline()
	aas = f.readline().split()
	for line in f:
		split = line.split()
		res = split.pop(0)
		BL50[res] = {}
		for i in range(len(split)):
			BL50[res][aas[i]] = int(split[i])
	BL50['-'] = {}
	for aa in aas:	#use -5 gap penalty (same as least conserved residues)
		BL50['-'][aa] = -5
		BL50[aa]['-'] = -5
 

def get_cons(spec, phos = 'default', nres=13, seed = 2122150, scoredict = BL50, Print = True):

	np.random.seed(seed)

	with open('blast/H.sapiens_'+spec+'.out', 'r') as f:
		seqh = []
		seqt = []
		indsh, indst = [], []
		for line in f:
			if line[:6] == 'ENSP00':
				seqh = seqh+[aa for aa in line.split()[1]]
				f.readline()
				seqt = seqt+[aa for aa in f.readline().split()[1] ]	

		n = 1
		for aa in seqh:
			if aa == '-': indsh.append('NA')
			else:
				indsh.append(n)
				n += 1

	if phos == 'random':
		phos = np.random.choice(range(1,915), size = nres, replace = False)

	#else: phos = [230, 249, 252, 356, 373, 608, 612, 780, 788, 795, 807, 811, 821, 826]
	else: phos = [249, 252, 356, 373, 608, 612, 780, 788, 795, 807, 811, 821, 826]

	phosh, phost = [seqh[indsh.index(p)] for p in phos], [seqt[indsh.index(p)] for p in phos]

	#print( [seqh[indsh.index(p)] for p in phos] )
	if Print: print('{:<15}'.format(spec.split('_')[-1])+'  '.join([seqt[indsh.index(p)] for p in phos]) )
	N = 0
	score = 0
	max_score = 0
	for i in range(len(phosh)):
		if phosh[i] == phost[i]: N += 1
		score += scoredict[phosh[i]][phost[i]]
		max_score += scoredict[phosh[i]][phosh[i]]
	return N, score, max_score #number of identical residues and similarity score
	
def compare_seqs(phos = 'default', add_para = False, nres=13, seed=2122150, Print =True):
	specs = ['C.jacchus', 'D.rerio', 'R.norvegicus', 'M.musculus', 'G.gallus', 'F.catus', 'M.mulatta', 'I.punctatus', 'M.gallopavo']
	sims = [97.5, 53, 86.5, 90.6, 72.4, 93.1, 96.8, 53.9, 71.2]

	if add_para:
		specs = specs+['RB1_H.sapiens_RBL1', 'RB1_H.sapiens_RBL2']
		sims = sims+[10, 5]
	specs = np.array(specs)
	sims = np.array(sims)
	args = np.argsort(-sims)
	sims = sims[args]/100
	specs = specs[args]

	p = 1.0
	for sim in sims: p = p*sim
	#print('p:', p)


	p = 1.0
	non_sim = 0
	score = 0
	max_score = 0
	for i, spec in enumerate(specs):
		N, new_score, new_max = get_cons(spec, phos=phos, nres=nres, seed=seed, Print = Print)
		new_p = 0.0
		for n in range(N, nres+1): new_p = new_p + comb(nres, N)*sims[i]**N*(1.0-sims[i])**(nres-N) #cumulative binomial distribution
	#print(new_p)
		p = p*new_p
		score += new_score
		max_score += new_max
	if Print: print('p:', p, 'score:', score, 'percent_max:', score/max_score*100)


	return p, score, score/max_score*100

def run_comparisons(N):
	seeds = range(N)
	ps = []
	scores = []
	for s in seeds:
		p, score, pmax = compare_seqs(phos = 'random', seed = s, Print=False)
		ps.append(-np.log10(p))
		scores.append(pmax)

	p, score, pmax = compare_seqs()

	reals = [-np.log10(p), pmax]
	xlabs = ['-log10(p)', 'relative similarity score']
	names = ['p', 'sim']
	data = [ps, scores]
	factors = [4, 2.5]
	nbins = int(N/1000)
	for i in range(2):
		plt.figure(figsize = (5,3.5))
		bins = np.linspace(min(data[i]), max(data[i]), nbins)
		plt.hist(data[i], bins = bins)
		plt.plot([reals[i], reals[i]], [0, int(N/nbins*factors[i])], 'k-')
		plt.xlabel(xlabs[i])
		plt.ylabel('frequency')
		plt.savefig('/local/data/public/ktj21/GIA3/figures/conservation_sim_'+names[i]+'.png',\
			dpi = 360, bbox_inches = 'tight')
		plt.show()

		print(min(data[i]), max(data[i]), np.mean(data[i]), reals[i])
		print('percentile:', len(np.array(data[i])[ np.array(data[i]) > reals[i] ]) / len(data[i]) )
	
if __name__ == '__main__':

	if len(sys.argv) > 1:
		N = int(sys.argv[1]) #can pick random residues
		run_comparisons(N)
	else:
		compare_seqs(add_para = True)

	

