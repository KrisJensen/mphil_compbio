from subprocess import call
from subprocess import Popen
import subprocess

specs = ['H.sapiens', 'M.musculus', 'M.mulatta', 'D.rerio']
isos = ['RB1', 'RBL1', 'RBL2']

for spec1 in specs:
	for iso1 in isos:
		for spec2 in specs:
			for iso2 in isos:
				if spec1+iso1 != spec2+iso2:
					#for each pair of non-identical isoforms, run ggsearch
					with open('blast/'+spec1+'_'+iso1+'_'+spec2+'_'+iso2+'.out', 'w') as f:
						call( ['/local/data/public/ktj21/programs/fasta36/bin/ggsearch36',\
							'-f', '-10', '-g', '-2',\
                	                        	spec1+'_'+iso1+'.pep.fa', '-subject', spec2+'_'+iso2+'.pep.fa'], stdout = f )




