#normalize mmseq data using deseq for different speices

tab = readRDS('tab_spec_isoforms.rds')

#genes for human and mouse
paralogs = rep(c('RB1', 'RBL1', 'RBL2'), 2)
names(paralogs) = c('ENST00000267163.5', 'ENST00000373664.7', 'ENST00000262133.10',
		'ENSMUST00000022701.6', 'ENSMUST00000029170.7', 'ENSMUST00000034091.7')
para_lens = c(4793, 4906, 5684, 4656, 4862, 4925)

tab = tab[names(paralogs),]

tab = tab*1000/para_lens #convert to RPK

#equate human and mouse genes
tab[1:3, c('Mmus_1', 'Mmus_2')] = tab[4:6, c('Mmus_1', 'Mmus_2')]
tab = tab[1:3,] #RPKM+1 and standardize genes
rownames(tab) = paralogs[1:3]
saveRDS(tab, 'tab_spec_parsed.rds')

tab = log2(tab+1)
par(mar = c(10,4,10,3))

cols = c('black', 'black', 'black', 'red', 'red', 'green', 'green')
pty = c(0, 1, 2, 4, 4, 8, 8)

#plot 2d graph of results
jpeg('species.jpg')
par(mar = c(10,4,10,3))
matplot(tab, type = 'o', lty = 2, col = cols, pch = pty, ylab = 'log2[RPKM+1]', axes = F)
legend(2.5,21, 'legend' = colnames(tab), col = cols, pch = pty, bty = 'n')
#legend(1.2,7.2, 'legend' = colnames(rep), col = cols, bty = 'n', pch = pty)
axis(1, at = seq(1, dim(tab)[1]), labels = rownames(tab), las=2)
axis(2)
dev.off()

#transform data for barplot; only use liver data
tabrep = tab[, 1:6]
tabrep[,1] = apply(tab[,4:5], 1, mean)
tabrep[,4] = apply(tab[,4:5], 1, sd)
tabrep[,2] = apply(tab[,1:3], 1, mean)
tabrep[,5] = apply(tab[,1:3], 1, sd)
tabrep[,3] = apply(tab[,6:7], 1, mean)
tabrep[,6] = apply(tab[,6:7], 1, sd)

mydata = data.frame(mean = as.numeric(t(tabrep[,1:3])), sd = as.numeric(t(tabrep[,4:6])))
mydata$names = c('Hep RB1', 'Hsap RB1', 'Mmus RB1',
		'Hep RBL1', 'Hsap RBL1', 'Mmus RBL1',
		'Hep RBL2', 'Hsap RBL2', 'Mmus RBL2')

jpeg('../../figures/species_bars.jpg')
lsize = 1.25
par(mar = c(10,4,10,3), cex.lab = lsize, cex.axis = lsize)
#plot means of expression
barCenters <- barplot(height = mydata$mean,
                  names.arg = mydata$names,
                  beside = true, las = 2,
                  ylim = c(0, max(mydata$mean+mydata$sd)),
		  col = rep(c('blue','green','red'), each=3),
                  cex.names = 0.75, xaxt = "n",
                  ylab = "log2[RPK+1]",
                  border = "black", axes = TRUE)

axis(1, at=barCenters, label=rep('', length(barCenters), lwd=2))
text(x = barCenters, y = par("usr")[3] - 1, srt = 45,
     adj = 1, labels = mydata$names, xpd = TRUE, cex = lsize)

#plot error bars
segments(barCenters, mydata$mean - mydata$sd, barCenters,
         mydata$mean + mydata$sd, lwd = 2)

arrows(barCenters, mydata$mean - mydata$sd, barCenters,
       mydata$mean + mydata$sd, lwd = 2, angle = 90,
       code = 3, length = 0.05)
dev.off()
