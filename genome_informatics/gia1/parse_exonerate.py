'''
Take output from exonerate and turn it into something useful
Provide keyword BCc/coli to align this species to the reciprocal species.
Use: python3 parse_exonerate.py BCc
'''

import sys
import pickle

def main(align = 'BCc'):

	if align == 'BCc': files, ind = ['BCc_genedict.pickled', 'exonerate_BCc_coli.out'], [4,5,1,0]
	elif align == 'coli': files, ind = ['coli_genedict.pickled', 'exonerate_coli_BCc.out'], [1,0,4,5]
	elif align == 'BCc_thresh': files, ind = ['BCc_genedict.pickled', 'exonerate_BCc_coli_thresh14.out'], [4,5,1,0]
	elif align == 'coli_thresh': files, ind = ['coli_genedict.pickled', 'exonerate_coli_BCc_thresh14.out'], [1,0,4,5]
	else: files, ind = align

	aligndict = {}

	genedict = pickle.load(open(files[0], 'rb'))

	for key in genedict.keys():
		aligndict[key] = []

	with open( files[1], 'r') as f:

		for line in f:
			split = line.split()
			if len(split) > 0 and 'Query:' == split[0]:
				q = split[ind[0]][ind[1]:]
				target = f.readline().split()
				if 'Target:' == target[0]:
					t = target[ind[2]][ind[3]:]
					f.readline()
					score = f.readline().split()[2]
					aligndict[q].append( [t, int(score)] )

	tot, aligned, multiple, unaligned, alignments, alignees = 0, 0, 0, 0, 0, []
	for key, item in aligndict.items():
		tot += 1
		if len(item) == 0: unaligned += 1
		else: aligned += 1
		if len(item) > 1: multiple += 1
		alignments += len(item)
		for i in item:
			if not i[0] in alignees: alignees.append(i[0])

	res = {'tot':tot, 'aligned':aligned, 'multiple':multiple, 'unaligned':unaligned, 'alignments':alignments, 'alignees':len(alignees)}

	print(res)

	return aligndict, res, alignees

if __name__ == '__main__':
	main( align = sys.argv[1] )


