import time
import subprocess


ts = []

t = time.time()
subprocess.call(['PyClone run_analysis_pipeline '+\
		'--in_files data/inputpyclone_noCN.tsv '+\
		'--working_dir pyclone_analysis_noCN '+\
		'--tumour_contents 0.8333333333333333333'],\
		shell=True)
t = time.time()-t
ts.append(t)
print 'time elsapsed:', t

t = time.time()
subprocess.call(['PyClone run_analysis_pipeline '+\
                '--in_files data/inputpyclone_0_7CN.tsv '+\
                '--working_dir pyclone_analysis_0_7CN '+\
		'--tumour_contents 0.8333333333333333333'],\
                shell=True)
t = time.time()-t
ts.append(t)
print 'time elsapsed:', t

print 'times:', ts


#PyClone run_analysis_pipeline --in_files data/PyClone_sim.tsv --working_dir pyclone_analysis_sim --tumour_contents 0.80
#PyClone build_mutations_file data/inputpyclone_noCN.tsv yaml/inputpyclone_noCN.yaml


