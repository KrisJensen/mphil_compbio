from subprocess import call

names = ['liver_emb_1', 'liver_emb_2',\
	#'liver_6_1', 'liver_6_2',\
	'liver_53_1', 'liver_53_2',\
	'breast-epithelium_51_1', 'breast-epithelium_51_2',\
	'spleen_51_1', 'spleen_51_2',\
	'colon_51_1', 'colon_51_2',\
	'adipose_53_1', 'adipose_53_2',\
	'pancreas_51_1', 'pancreas_51_2',\
	'heart_51_1', 'heart_51_2'\
	]
files = ['ENCFF001RNR', 'ENCFF001RNQ',\
	#'ENCFF001RUA', 'ENCFF001RTZ',\
	'ENCFF381FTA', 'ENCFF151LEE',\
	'ENCFF421LZH', 'ENCFF613GRK',\
	'ENCFF455UMQ', 'ENCFF361OLP',\
	'ENCFF097YEA', 'ENCFF593EMI',\
	'ENCFF460GEN', 'ENCFF222EUQ',\
	'ENCFF517RES', 'ENCFF285HWP',\
	'ENCFF855TSM', 'ENCFF971MYL'\
	]

for i in range(len(names)):
	print('new file', names[i])
	call( ['wget', '-O', names[i]+'.fastq.gz',\
		'https://www.encodeproject.org/files/'+files[i]+'/@@download/'+files[i]+'.fastq.gz'\
		] )

#wget -O liver_6_1 https://www.encodeproject.org/files/ENCFF001RUA/@@download/ENCFF001RUA.fastq.gz
