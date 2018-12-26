'''
After running exonerate, we can analyze the output file
to find aligned genes and then find the corresponding GO terms
'''

import sys
from parse_exonerate import main
import pickle

gene_GO = pickle.load(open('geneGOs.pickled', 'rb'))
GOfunctions = pickle.load(open('GOFunctions.pickled', 'rb'))

aligndict, res, alignees = main(align = sys.argv[1])

counts = {}

for bcc, colis in aligndict.items():
	print(bcc, colis)
	for coli in colis: #find each aligned E.coli gene
		#print(bcc, colis, coli)
		try:
			for go in gene_GO[coli[0]]: #get GO term for each aligned E.coli gene
				if not go in counts.keys(): counts[ go ] = 1 
				else: counts[ go ] += 1 #add count
		except KeyError:
			print(coli[0], 'not found')

ordered = list(sorted(counts, key = counts.get, reverse = True)) #sort GO terms by counts
vals = [counts[key] for key in ordered]

for i in range(8): print(vals[i], ordered[i], GOfunctions[ordered[i]], '\n') 
	



