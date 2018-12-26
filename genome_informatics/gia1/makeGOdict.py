'''parses the go-basis.obo file an dumps a pickled dict of GO terms and functions'''

import pickle

terms = {}

with open('go-basic.obo', 'r') as f:

	for line in f:
		if line[:3] == 'id:':
			go = line.split()[1]
			newline = f.readline()
			if 'name:' == newline[:5]:
				terms[go] = ' '.join(newline.split()[1:])

gene_functions = {}
gene_GO = {}
with open('gene_association.ecocyc', 'r') as f:
	for i in range(19): f.readline()
	for line in f:
		split = line.split()
		genes = [ split[2] ]
		go = split[3] 

		for word in split:
			if '|' in word:
				genes = genes+word.split('|')

		if go in terms.keys():
			for gene in genes:
				if gene in gene_functions.keys() and not go in gene_GO[gene]:
					if go[:2] == 'GO':
						gene_GO[gene].append( go )
						gene_functions[gene].append( terms[go] )
				else:
					if go[:2] == 'GO':
						gene_GO[gene] = [ go ]
						gene_functions[gene] = [ terms[go] ]


pickle.dump(terms, open('GOFunctions.pickled', 'wb'))
pickle.dump(gene_functions, open('geneFunctions.pickled', 'wb'))			
pickle.dump(gene_GO, open('geneGOs.pickled', 'wb'))

