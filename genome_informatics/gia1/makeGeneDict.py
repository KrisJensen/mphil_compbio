'''Construct simple dictionary of genes from B. a'''

import pickle

genes = {}

with open('Buchnera_aphidicola_bcc.ASM9096v1.cdna.all.fa', 'r') as f:

	for line in f:
		if '>' == line[0]:
			split = line.split()
			genes[split[3][5:]] = split[6][12:]

pickle.dump(genes, open('BCc_genedict.pickled', 'wb'))


genes = {}

with open('Escherichia_coli_str_k_12_substr_mds42.GCA_000350185.1.23.cdna.all.names.fa', 'r') as f:

        for line in f:
                if '>' == line[0]:
                        split = line.split()
                        genes[split[0][1:]] = ''.join(split[6:])

pickle.dump(genes, open('coli_genedict.pickled', 'wb'))


