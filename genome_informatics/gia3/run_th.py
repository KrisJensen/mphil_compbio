'''
given a list of sample basenames, aligns these
to a reference chromosome using tophat
'''

from subprocess import call
import sys

names = []
for arg in sys.argv[1:]: #construct list of samples
	names.append( (arg+'_1_val_1', arg+'_2_val_2') )

for files in names:
	name = files[0][:-8]
	print('new sample', name)
	if not 'Mmus' in name: #align to human chromosome 13
		call( ['tophat', '--max-multihits', '1', '-p', '12', '-o', '../bam/dna/'+name+'_th',\
		'dna/GRCh38/Homo_sapiens.GRCh38.dna.chromosome.13',\
		'fastq/'+files[0]+'.fq', 'fastq/'+files[1]+'.fq'\
		] )

	else: #align to mouse chromosome 14
		call( ['tophat', '--max-multihits', '1', '-p', '12', '-o', '../bam/dna/'+name+'_th',\
		'dna/GRCm38/Mus_musculus.GRCm38.dna.chromosome.14',\
		'fastq/'+files[0]+'.fq', 'fastq/'+files[1]+'.fq'\
		] )

