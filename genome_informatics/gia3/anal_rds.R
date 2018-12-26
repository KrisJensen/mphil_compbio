#plot some pretty picture from gene-level analysis

require(DESeq2)

dds = readRDS("dds_genes.rds") #load dataset

#genes to look at (all on chr13)
ref_genes = c('RB1', 'BRCA2', 'SOX21', 'HTR2A', 'CCNA1')

names(ref_genes) = c('ENSG00000139687',
		'ENSG00000139618',
		'ENSG00000125285',
                'ENSG00000102468',
                'ENSG00000133101')

reps = c('hepatocyte_1', 'hepatocyte_2',
        'epithelium_1', 'epithelium_2',
        'cardiac_1', 'cardiac_2',
        'neural-progenitor_1', 'neural-progenitor_2',
	'liver_emb', 'liver_6', 'liver_53',
	'adipose_53', 'spleen_51', 'pancreas_51',
	'colon_51', 'heart_51', 'breast-epithelium_51')

ref_lens = c(4793, 11986, 5429, 2778, 1879) #just pick consensus transcripts

dds_ref = dds[names(ref_genes), reps]

#pick out data of interest
RB1 = counts(dds, normalized=TRUE)['ENSG00000139687',]
ref = counts(dds_ref, normalized=TRUE)
rownames(ref) = ref_genes
print(ref)

#write some information to files
write.table(sizeFactors(dds), file = 'size_factors.txt')
write.table(RB1, file = 'geneCounts.txt')
write.table(ref, file = 'geneCounts_all.txt')

#normalize and log-transform
ref = 1000*ref/ref_lens
rep = log2(ref+1)

means = apply(rep, 1, mean); sds = apply(rep, 1, sd)
print(means)
print(sds)

#jpeg('genes.jpg')
par(mar = c(10,4,10,3))

#plot 2d graph of expression
linerep = rep[,1:8]
cols = c('black', 'black', 'red', 'red', 'green', 'green', 'blue', 'blue')
pty = c(0,0, 4, 4, 8, 8, 1, 1)
matplot(linerep, type = 'o', lty = 2, col = cols, pch = pty, ylab = 'log2[RPK+1]', axes = F)
legend(1,7.5, 'legend' = colnames(linerep), col = cols, pch=pty, bty = 'n')
axis(1, at = seq(1, dim(rep)[1]), labels = rownames(linerep), las=2)
axis(2)
#dev.off()

r = cor(rep[1,], rep[5,]) #correlation between RB1 and CCNA1
cat('correlation coefficient', r, '\n')

cat('standard deviations', sd(rep[1,]), sd(rep[5,]), '\n')


#make barplot of RB1 or CCNA1 expression
prot = 'CCNA1'
name = 'RB1_genes_bars'
ind = 1

if (prot == 'CCNA1'){
	ind = 5
	name = 'CCNA1_genes_bars'
}

cat('analyzing ', prot, '\n')

tabrep = rep[1,1:10]

#reformat data
for (i in 1:4){
	cat('new i ', i, '\n')
        tabrep[i] = mean(rep[ind, (2*i-1):(2*i) ])
        tabrep[i+5] = sd(rep[ind, (2*i-1):(2*i) ])
}

cat('summarized data\n')

tabrep[5] = mean(rep[ind, 9:17])
tabrep[10] = sd(rep[ind, 9:17])

cat('finished organs\n')

mydata = data.frame(mean = as.numeric(t(tabrep[1:5])), sd = as.numeric(t(tabrep[6:10])))
mydata$names = c('hepatocyte', 'epithelium', 'cardiac', 'neural', 'organs')
lsize = 1.25
jpeg(paste0('../../figures/', name, '.jpg'))
par(mar = c(10,4,10,3), cex.axis = lsize, cex.lab = lsize)
#plot means
barCenters <- barplot(height = mydata$mean,
                  names.arg = mydata$names,
                  beside = true, las = 2,
                  ylim = c(0, max(mydata$mean+mydata$sd)+0.2),
                  col = c('blue','green','red','yellow', 'pink'),
                  cex.names = 0.75, xaxt = "n",
                  ylab = "log2[RPK+1]",
                  border = "black", axes = TRUE)
		  #main = paste(prot, 'expression'))

axis(1, at=barCenters, label=rep('', length(barCenters), lwd=2))
text(x = barCenters, y = par("usr")[3] - 1, srt = 45,
     adj = 1, labels = mydata$names, xpd = TRUE, cex=lsize)
#plot error bars
segments(barCenters, mydata$mean - mydata$sd, barCenters,
         mydata$mean + mydata$sd, lwd = 2)

arrows(barCenters, mydata$mean - mydata$sd, barCenters,
       mydata$mean + mydata$sd, lwd = 2, angle = 90,
       code = 3, length = 0.05)
dev.off()


