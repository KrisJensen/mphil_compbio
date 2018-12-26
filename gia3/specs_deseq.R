#perform differential expression analysis for mouse vs human

require(DESeq2)

#load counts
counts = readRDS("tab_spec_parsed.rds")
samples = list(c('hepatocyte_1','hepatocyte_2'),
		c('liver_emb', 'liver_6', 'liver_53'),
		c('Mmus_1', 'Mmus_2'))
types = c('Hep', 'Hsap', 'Mmus')

paralogs = c('RB1', 'RBL1', 'RBL2')
genes = c('ENST00000267163.5', 'ENST00000373664.7', 'ENST00000262133.10')

for (k in 1:3){ #for each paralog

pvals = matrix(, 3, 3) #result matrix
rownames(pvals)=colnames(pvals)=types

for (i in 1:3){
	for (j in 1:3){ #analzye each pair of samples
		if (i != j){
		sampnames = c(samples[[i]], samples[[j]])
		lens = c(length(samples[[i]]), length(samples[[j]]))

		newcounts = counts[,sampnames] #pick out relevant data
		colData = matrix( c( rep(types[[i]], lens[1]),
					rep(types[[j]], lens[2]) ), sum(lens), 1 )
		rownames(colData) = sampnames
		colnames(colData) = 'Tissue'
		#construct DEseq dataset
		dds <- DESeqDataSetFromMatrix(round(newcounts),
			colData = colData, design = ~ Tissue)
		normalizationFactors(dds) = matrix(1, 3, sum(lens)) #already normalized
		dds <- estimateDispersions( dds )
		dds = nbinomWaldTest(dds) #conduct Wald's test
		res = results(dds)
		p = res[paralogs[k], 'pvalue']
		pvals[i,j] = p #store result
		}
	}
}

pvals = pvals * 3 * 3 #bonferroni correction
cat('\n\n\n', paralogs[k], '\n')
print(pvals) #print result
}

