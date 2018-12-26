'''
runs exonerate to align B. aphidicola genes to E. coli genes for a range of
different percent thresholds.
'''

import numpy as np
from subprocess import call

for i in np.arange(0.5, 100.5, 0.5):
	with open('exonerate_screen_2/exonerate_BCc_coli_p'+str(i)+'.out', 'w') as f:
		call( [ 'exonerate', '--query', 'Buchnera_aphidicola_bcc.ASM9096v1.cdna.all.fa',\
			'--target', 'Escherichia_coli_str_k_12_substr_mds42.GCA_000350185.1.23.cdna.all.names.fa',\
			'--percent', str(i), '--score', '20'], stdout = f )
