'''
Given a list of sample basenames, processes the bam file generated
by align_mm.py using mmseq to map the alignments to genes
creates .mmseq files with reads per transcripts
'''

from subprocess import call
import sys

names = []
for arg in sys.argv[1:]: #construct list of sample names
	names.append( (arg+'_1', arg+'_2') )

for files in names:
	name = files[0][:-2]
	print('new sample', name)
	#run bam2hits from mmseq
	with open('../bam/cdna/'+name+'_mm/'+name+'.hits', 'w') as f:
		if not 'Mmus' in name: #map to human genes
			call( ['bam2hits', 'cdna/GRCh38/Homo_sapiens.GRCh38.cdna.all.fa', '../bam/cdna/'+name+'_mm/'+name+'.namesorted.bam'],\
			stdout = f )
		else: #map to mouse genes
			call( ['bam2hits', 'cdna/GRCm38/Mus_musculus.GRCm38.cdna.all.fa', '../bam/cdna/'+name+'_mm/'+name+'.namesorted.bam'],\
			stdout = f )

	#obtain expression estimates
	call( ['mmseq', '../bam/cdna/'+name+'_mm/'+name+'.hits', '../bam/cdna/'+name+'_mm/'+name] )

