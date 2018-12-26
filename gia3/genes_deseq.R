#perform gene-level differential expression analysis

require(DESeq2)

counts = readRDS("matrix_genes.rds")

#only consider data where we have duplicates
types = c('hepatocyte', 'epithelium', 'cardiac', 'neural-progenitor')

pvals = matrix(, 4, 4) #initialize result matrix
rownames(pvals)=colnames(pvals)=types

for (i in 1:4){ 
	for (j in 1:4){ #do all pairwise comparisons
		if (i != j){
		type1 = types[i]; type2 = types[j]
		sampnames = c(paste0(type1, '_1'), paste0(type1,'_2'),
                                paste0(type2, '_1'), paste0(type2,'_2'))

		newcounts = counts[,sampnames] #get counts for these samples
		colData = matrix( c(type1, type1, type2, type2), 4, 1)
		rownames(colData) = sampnames
		colnames(colData) = 'Tissue'
		#construct deseq dataset
		dds <- DESeqDataSetFromMatrix(newcounts,
			colData = colData, design = ~ Tissue)
		dds = DESeq(dds) #analyze results
		res = results(dds)
		p = res['ENSG00000139687', 'pvalue']
		pvals[i,j] = p #store result
		newcounts = counts(dds, normalize=TRUE)
                if (i==2){
                        cat(newcounts['ENSG00000139687',])
                        cat('epithelium - ',types[j],': ',
                        (sum(newcounts['ENSG00000139687',1:2])-sum(newcounts['ENSG00000139687',3:4]))/sum(newcounts['ENSG00000139687',1:2]),
                        '\n\n')
                }

		}
	}
}

pvals = pvals * 6 #bonferroni correction
print(pvals)
