from subprocess import call


names = ['hepatocyte_1_1', 'hepatocyte_1_2',\
	'hepatocyte_2_1', 'hepatocyte_2_2',\
	'epithelium_1_1', 'epithelium_1_2',\
	'epithelium_2_1', 'epithelium_2_2',\
	'cardiac_1_1', 'cardiac_1_2',\
	'cardiac_2_1', 'cardiac_2_2',\
	'neural-progenitor_1_1', 'neural-progenitor_1_2',\
	'neural-progenitor_2_1', 'neural-progenitor_2_2',\
	'Mmus_1_1', 'Mmus_1_2',\
	'Mmus_2_1', 'Mmus_2_2',\
	]
files = ['ENCFF245VTB','ENCFF369QXD','ENCFF653EZE','ENCFF644RGE',\
	'ENCFF464PFE', 'ENCFF987MZF','ENCFF846YCY','ENCFF751KJR',\
	'ENCFF353SPR', 'ENCFF197MZY','ENCFF030WXU','ENCFF724FYF',\
	'ENCFF939FVE', 'ENCFF201WLO', 'ENCFF996KMK', 'ENCFF726PAY',\
	'ENCFF001RTN', 'ENCFF001RTM', 'ENCFF001RTD', 'ENCFF001RTA',\
	]

for i in range(len(names)):
	print('new file', names[i])
	call( ['wget', '-O', names[i]+'.fastq.gz',\
		'https://www.encodeproject.org/files/'+files[i]+'/@@download/'+files[i]+'.fastq.gz'\
		] )

#wget -O liver_6_1 https://www.encodeproject.org/files/ENCFF001RUA/@@download/ENCFF001RUA.fastq.gz
