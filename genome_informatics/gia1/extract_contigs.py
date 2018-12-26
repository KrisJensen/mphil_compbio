'''Given a directory name dir and an integer n, finds the n largest contigs in dir and puts them in separate fasta files
for subsequent BLAST'''

import sys
import numpy as np

directory, n = sys.argv[1], int(sys.argv[2])

lengths = []

with open(directory+'/stats.txt', 'r') as f:
	f.readline() #read header
	for line in f: lengths.append( int(line.split()[1]) ) #store lengths

lengths = list(reversed(np.sort(lengths)))[:n]

print('finding contigs of lengths', lengths)

with open(directory+'/contigs.fa', 'r') as f:

	for line in f:
		for l in lengths:
			if 'length_'+str(l) in line:
				print('found length', l)
				with open(directory+'/contig'+str(lengths.index(l))+'.fa', 'w') as fout:
					finished = False
					L = 61
					while L == 61: #standard length of line in velvet .fa file. Last line is less than this.
						newline = f.readline()
						L = len(newline)
						#print(newline, L)
						if not 'NODE' in newline: fout.write(newline) #if by chance N%60 = 0, this stops us writing
					


