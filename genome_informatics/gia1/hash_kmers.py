'''Hashes our sequencing data for a number of different kmer lengths for further processing with velvetg'''

import sys
from subprocess import call


for i in range(17,99,2): #hardcode kmer lengths. Below 17, we don't assemble a sensible graph and 97 is almost our full read length
	call(['/local/data/public/ktj21/programs/velvet/velveth', 'kmer_'+str(i), str(i), '-fastq', '-separate',\
	'/local/data/public/genome_informatics_2018/assignments/assignment_1/read1.fq',\
	'/local/data/public/genome_informatics_2018/assignments/assignment_1/read2.fq']) #run velveth






