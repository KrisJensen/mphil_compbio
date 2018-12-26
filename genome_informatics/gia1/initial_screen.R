#reproduce the analysis of coverage from the velvet manual


library(plotrix)

data = read.table("./dist_kmer_cut7/stats.txt", header=TRUE)
cov = data$short1_cov
lgth = data$lgth

pdf(file = './figs/hist_kmer57_cut7.pdf')

hist(cov, xlim=range(0,50), breaks=0:10000, #plot raw histogram
     main = 'Unweighted Coverage cut_cutoff 7', xlab = 'Coverage')

dev.off()


pdf(file = './figs/weighthist_kmer57_cut7.pdf')

weighted.hist(cov, lgth, breaks=0:50, #plot weighted histogram
      main = 'Weighted Coverage cut_cutoff 7', xlab = 'Coverage')

dev.off()

#weighted.hist(data$short1_cov, data$lgth, breaks=0:50)
