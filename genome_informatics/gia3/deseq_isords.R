
require(DESeq2)

#dds = readRDS("dds_temp.rds")
dds = readRDS("dds_isoforms.rds")

#rld <- rlog(dds, blind=FALSE)
#plotPCA(rld, intgroup='Tissue')


ref_genes = c('RB1-201', 'RB1-209', 'RB1-208', 'RB1-207', 'RB1-202', 'RB1-206', 'RB1-203', 'RB1-205', 'RB1-204')
ref_lens = c(4793, 4629, 633, 329, 480, 658, 617, 1490, 1097)

names(ref_genes) = c('ENST00000267163.5',
		'ENST00000650461.1',
                'ENST00000646097.1',
                'ENST00000643064.1',
                'ENST00000467505.5',
		'ENST00000531171.5',
		'ENST00000480491.1',
		'ENST00000525036.1',
		'ENST00000484879.1')

paralogs = c('RB1-201', 'RBL1-202', 'RBL2-201')
names(paralogs) = c('ENST00000267163.5', 'ENST00000373664.7', 'ENST00000262133.10')
para_lens = c(4793, 5684, 4906)

reps = c('hepatocyte_1', 'hepatocyte_2',
	'epithelium_1', 'epithelium_2',
	'cardiac_1', 'cardiac_2',
	'neural-progenitor_1', 'neural-progenitor_2')
specs = c('hepatocyte_1', 'hepatocyte_2',
	'Mmus_1', 'Mmus_2',
	'liver_53')

#dds_spec = estimateDispersions(dds[names(paralogs), specs])
dds_rep = estimateDispersions(dds[names(ref_genes), reps])

#res_spec = nbinomTest(dds_spec, "hepatocye", "Mmus")
res_rep = nbinomTest(dds_rep, "hepatocyte","epithelium")


rep = counts(dds_rep, normalized=TRUE)
rownames(rep) = ref_genes

#spec = counts(dds_spec, normalized=TRUE)
#rownames(spec) = paralogs

print(rep)
write.table(sizeFactors(dds), file = 'size_factors_transcripts_reps.txt')
write.table(rep, file = 'transcriptCounts_reps.txt')

rep = log2(rep+1)

means = apply(rep, 1, mean); sds = apply(rep, 1, sd)
print(means)
print(sds)

cols = c('black', 'blue', 'red', 'green', 'yellow')
matplot(rep, type = 'l', lty = 2, col = cols, ylab = 'log2[reads+1]', axes = F)
legend(c(5,5), 'legend' = colnames(ref), col = cols, fill = cols, bty = 'n')
axis(1, at = seq(1, dim(ref)[2]), labels = rownames(rep))

#matplot(refp, type = 'l', lty = 2, col = cols, ylab = 'log2[reads+1]', axes = F)


