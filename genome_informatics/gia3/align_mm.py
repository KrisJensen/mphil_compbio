'''
given a list of fastq basenames as arguments, aligns
the corresponding fastq files to the relevant reference cdna
for further processing using mmseq
'''

from subprocess import call
import sys

names = []
for arg in sys.argv[1:]: #construct list of fq filenames
	names.append( (arg+'_1_val_1', arg+'_2_val_2') )

for files in names:
	name = files[0][:-8]
	print('new sample', name)
	if not 'Mmus' in name: #align to human cdna
		call( ['bowtie -a --best --strata -S -m 100 -X 500 --chunkmbs '+\
		'256 -p 12 cdna/GRCh38/Homo_sapiens.GRCh38.cdna.all -1 fastq/'+files[0]+'.fq -2 fastq/'+files[1]+'.fq '+\
		'| samtools view -F 0xC -bS - | samtools sort -n - '+\
		'../bam/cdna/'+name+'_mm/'+name+'.namesorted'\
		], shell=True )

	else: #align to mouse cdna
		call( ['bowtie -a --best --strata -S -m 100 -X 500 --chunkmbs '+\
		'256 -p 12 cdna/GRCm38/Mus_musculus.GRCm38.cdna.all -1 fastq/'+files[0]+'.fq -2 fastq/'+files[1]+'.fq '+\
		'| samtools view -F 0xC -bS - | samtools sort -n - '+\
		'../bam/cdna/'+name+'_mm/'+name+'.namesorted'\
		], shell=True )

