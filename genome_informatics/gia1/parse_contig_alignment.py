'''counts how many of the BCc genes we find in our assembly'''

def get_genes():
	genes = []

	with open('exonerate_BCc_contig.out', 'r') as f:

		for line in f:
			split = line.split()
			if len(split) > 0 and split[0] == 'Query:':
				if not split[4] in genes:
					genes.append(split[4])

	print('found', len(genes), 'genes at 95% threshold')

	return(genes)

if __name__ == '__main__':
	get_genes()

