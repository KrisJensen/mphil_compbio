'''
given a list of sample basenames,
trims the corresponding fastq files using trimgalore
'''

from subprocess import call
import sys

names = []
for arg in sys.argv[1:]: #construct list of samples
	names.append( (arg+'_1', arg+'_2') )


for files in names:
	name = files[0][:-2]
	print('new sample', name)
	#run trimgalore using parameters suggested by mmseq
	call( ['trim_galore', '-q', '15', '--stringency', '3', '-e', '0.05', '--length', '36',\
		'--trim1', '--paired', files[0]+'.fastq', files[1]+'.fastq'\
		] )

