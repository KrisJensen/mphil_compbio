from subprocess import call
from subprocess import Popen
import subprocess

#define species of interest for which we have data
specs = ['H.sapiens', 'C.elegans', 'D.rerio', 'R.norvegicus', 'M.musculus',\
	'G.gallus', 'F.catus', 'M.mulatta', 'C.jacchus', 'I.punctatus', 'M.gallopavo']


def run_hmmer(spec):
	'''given a species, runs a hmmer alignment to pfam'''
	call( ['/local/data/public/genome_informatics_2018/programs/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmscan', '--cpu', '60',\
                '--domtblout', 'hmmer/'+spec+'_tbl.txt', '/local/data/public/genome_informatics_2018/assignments/assignment_2/pfam_database/Pfam-A.hmm',\
                spec+'.pep.fa'] )

for spec in specs+['RBL1', 'RBL2']: #also align homologs RBL1 and RBL2
	run_hmmer(spec) #align the proteins

