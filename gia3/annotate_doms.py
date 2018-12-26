'''parses the hmmer output files to generate a list
with all proteins and domain annotations'''

import pickle
import numpy as np
import sys

#For each protein;
#if domain aligned to < thresh: add to list, annotate as domain, store GO

def get_best_doms(doms):
	'''given a list of domains passing threshold,
	picks out the most significant ones and assembles them without overlap'''
	GO_num = pickle.load(open('GO_nums.pickled', 'rb')) #get GO terms
	doms_out = []
	prot_GOs = []
	targets, evals, ranges, evals_num = [np.array([]) for i in range(4)]
	for dom in doms: #store E values
		evals_num = np.append(evals_num, float(dom['Evalue']))

	args = np.argsort(evals_num)
	doms = np.array(doms)[args] #sort by Evalue

	used = set() #keep track of residues aligned to
	for dom in doms:
		newrange = range(int(dom['range'][0]), int(dom['range'][1])+1)
		if len(used.intersection(newrange)) < 10: #only add domain if it doesn't overlap with already annotated domains
			used.update(newrange)
			t = dom['target']
			if t in GO_num.keys(): #add GO terms
				GOs = GO_num[t]
				dom['GOs'] = GOs
				for GO in GOs:
					if not GO in prot_GOs: prot_GOs.append(GO)
			else: dom['GOs']=[]
			doms_out.append(dom)

	return doms_out, prot_GOs

def parse_domtbl(spec, thresh = 10**(-1), i=0):
	'''parses a hmmer outputfile'''
	prots = {}
	q0 = 'NA'
	with open('hmmer/'+spec+'_tbl.txt', 'r') as f:
		for line in f:
			if not line[0] == '#':
	
				split = line.split() #each line is a new alignment
				t, q, q_a, E = split[1].split('.')[0], split[3], split[4], split[12]

				if q != q0: #if we're looking at a new query protein
					if not q0 == 'NA':
						best_doms, GOs = get_best_doms(doms) #parse previous query
						prots[q0]['doms'] = best_doms
						prots[q0]['GOs'] = GOs

					qlen = split[5]
					prots[q] = {'length':qlen}
					used = set()
					q0 = q
					doms = []


				if float(E) < thresh: #if same protein and alignment significant, add to list
					doms.append({'range': [split[17], split[18]], 'target':t, 'Evalue':E})

		#add last protein
		best_doms, GOs = get_best_doms(doms)
		prots[q]['doms'] = best_doms
		prots[q]['GOs'] = GOs

	if i==0: method = 'w'
	else: method = 'a'	
			
	with open('hmmer/doms_parsed', method) as f: #write results
		for q, align in prots.items():
			print(q, align['length'], str(len(align['doms'])), align['GOs'])
			f.write('>'+q+'   qlen:'+align['length']+'   doms:'+str(len(align['doms']))+'   GOs:'+','.join(align['GOs'])+'\n')
			
			for dom in align['doms']:
				f.write('   target:'+dom['target']+'   range:'+','.join(dom['range'])+'   Evalue:'+dom['Evalue']+'   GOs:'+','.join(dom['GOs'])+'\n')
			

	pickle.dump(prots, open('hmmer/'+spec+'_domains_m1.pickled', 'wb')) #save alignment

if __name__ == '__main__':

	specs = ['H.sapiens', 'C.jacchus', 'D.rerio', 'R.norvegicus', 'M.musculus',
	        'G.gallus', 'F.catus', 'M.mulatta', 'I.punctatus', 'M.gallopavo', 'RBL1', 'RBL2']
	sims = [100, 97.5, 53, 86.5, 90.6,
	        72.4, 93.1, 96.8, 53.9, 71.2, 10, 5]
	specs = np.array(specs)
	sims = np.array(sims)
	args = np.argsort(-sims)
	sims = sims[args]/100
	specs = specs[args] #sort by similarity

	for i, spec in enumerate(specs):
		parse_domtbl(spec, i=i) #parse output for each protein

