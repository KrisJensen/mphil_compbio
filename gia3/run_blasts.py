'''doesn't actually run blasts, runs ggsearch36 for global alignment'''

from subprocess import call
from subprocess import Popen
import subprocess

#define proteins to analyze
specs = ['H.sapiens', 'C.elegans', 'D.rerio', 'R.norvegicus', 'M.musculus', 'G.gallus', 'F.catus', 'M.mulatta', 'C.jacchus', 'I.punctatus', 'M.gallopavo']

for spec1 in specs:
	for spec2 in specs:
		if spec1 != spec2: #for each pair of different species, run ggsearch36
			with open('blast/'+spec1+'_'+spec2+'.out', 'w') as f:
				call( ['/local/data/public/ktj21/programs/fasta36/bin/ggsearch36',\
					'-f', '-10', '-g', '-2',\
                                        spec1+'.pep.fa', '-subject', spec2+'.pep.fa'], stdout = f )




