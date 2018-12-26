'''
runs exonerate to align B. aphidicola genes to E. coli genes for a range of
different raw score thresholds. This isn't the best way of doing it but just
repeated the analysis for percent score.
'''

import numpy as np
from subprocess import call

for i in np.arange(20, 201, 1):
	with open('exonerate_screen_raw/exonerate_BCc_coli_r'+str(i)+'.out', 'w') as f:
		call( [ 'exonerate', '--query', 'Buchnera_aphidicola_bcc.ASM9096v1.cdna.all.fa',\
			'--target', 'Escherichia_coli_str_k_12_substr_mds42.GCA_000350185.1.23.cdna.all.names.fa',\
			'--score', str(i)], stdout = f )
