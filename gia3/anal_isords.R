#plot mmseq output data
require(DESeq2)

dds = readRDS("dds_isoforms.rds")

#pick out RB1 transcripts
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

paralogs = c('RB1', 'RBL1', 'RBL2')
names(paralogs) = c('ENST00000267163.5', 'ENST00000373664.7', 'ENST00000262133.10')
para_lens = c(4793, 4906, 5684)

#we can pick out different groups of samples
reps = c('hepatocyte_1', 'hepatocyte_2',
	'epithelium_1', 'epithelium_2',
	'cardiac_1', 'cardiac_2',
	'neural-progenitor_1', 'neural-progenitor_2')
specs = c('hepatocyte_1', 'hepatocyte_2',
	'Mmus_1', 'Mmus_2',
	'liver_53')

organs = c('adipose_53', 'spleen_51', 'pancreas_51', 'colon_51', 'heart_51', 'breast-epithelium_51',
        'liver_emb', 'liver_6', 'liver_53')

dds_rep = dds[names(ref_genes), c(reps, organs)]
dds_para = dds[names(paralogs), reps]
dds_allpara = dds[names(paralogs), c(reps, organs)]

rep = counts(dds_rep, normalized=TRUE); rownames(rep) = ref_genes
allpara = counts(dds_allpara, normalized=TRUE); rownames(allpara) = paralogs

#make boxplot
lsize = 1.3
jpeg('../../figures/transcripts_boxplot.jpg')
par(mar = c(10,4,10,3), cex.lab = lsize, cex.axis = lsize)

boxplot.matrix(log2(t(rep+1)), ylab='log2[normalized reads + 1]', las=2)
dev.off()

fracs = apply(rep[1:2,],2, sum)/apply(rep,2, sum)
print(fracs)
print(mean(fracs))

rep = rep[,reps]

#normalize data
rep = rep*1000/ref_lens #reads per kilobase
rep = counts(dds_para, normalized=TRUE); rownames(rep) = paralogs; rep = rep*1000/para_lens 

allpara = log2(allpara*1000/para_lens+1)
for (i in 1:3){
	jpeg(paste0(paralogs[i], '_variation_bars.jpg'))
	barplot(allpara[i,], las=2)
	dev.off()
	cat('\nstats for ', paralogs[i], ' : ', mean(allpara[i,]), ' ', sd(allpara[i,]), '\n')
}

#spec = counts(dds_spec, normalized=TRUE)
#rownames(spec) = paralogs

print(rep)
write.table(rep, file = 'transcriptCounts_reps.txt')

rep = log2(rep+1) #log transform for plotting

means = apply(rep, 1, mean); sds = apply(rep, 1, sd)
print(means)
print(sds)

par(mar = c(10,4,10,3))

#make plots look nice and plot 2d graph
cols = c('black', 'black', 'red', 'red', 'green', 'green', 'blue', 'blue')
pty = c(0,0, 4, 4, 8, 8, 1, 1)

#pty = 1:9; cols = 'black'
#jpeg('transcripts.jpg')
jpeg('transcripts_paralogs.jpg')
par(mar = c(10,4,10,3))
matplot(rep, type = 'o', lty = 2, col = cols, pch = pty, ylab = 'log2[RPK+1]', axes = F)
#legend(2.8,9, 'legend' = colnames(rep), col = cols, pch = pty, bty = 'n')
legend(1.2,7.8, 'legend' = colnames(rep), col = cols, bty = 'n', pch = pty)
axis(1, at = seq(1, dim(rep)[1]), labels = rownames(rep), las=2)
axis(2)
dev.off()

par(mar = c(10,4,10,3))
barplot(log2(2^rep[1,]+2^rep[2,]), las = 2, ylab='log2[RPK+1]')


#matplot(refp, type = 'l', lty = 2, col = cols, ylab = 'log2[reads+1]', axes = F)

#transform data for bar plot
newrep = rep
newrep[1,] = log2(2^newrep[1,]+2^newrep[2,]) #consider two major transcripts
tabrep = newrep[1,]

for (i in 1:4){
	tabrep[i] = mean(newrep[1, (2*i-1):(2*i) ])
	tabrep[i+4] = sd(newrep[1, (2*i-1):(2*i) ])
}

mydata = data.frame(mean = as.numeric(t(tabrep[1:4])), sd = as.numeric(t(tabrep[5:8])))
mydata$names = c('hepatocyte', 'epithelium', 'cardiac', 'neural')

#jpeg('RB1_transcripts_bars.jpg')
par(mar = c(10,4,10,3))
#plot means
barCenters <- barplot(height = mydata$mean,
                  names.arg = mydata$names,
                  beside = true, las = 2,
                  ylim = c(0, max(mydata$mean+mydata$sd)),
                  col = c('blue','green','red','yellow'),
                  cex.names = 0.75, xaxt = "n",
                  ylab = "log2[RPK+1]",
                  border = "black", axes = TRUE)

text(x = barCenters, y = par("usr")[3] - 1, srt = 45,
     adj = 1, labels = mydata$names, xpd = TRUE)

#plot error bars
segments(barCenters, mydata$mean - mydata$sd, barCenters,
         mydata$mean + mydata$sd, lwd = 1.5)

arrows(barCenters, mydata$mean - mydata$sd, barCenters,
       mydata$mean + mydata$sd, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)
dev.off()
