#script for running a ccube calclation and plotting the result
require('ccube')
require('foreach')
require('dplyr') 
library('reticulate')

args = commandArgs(trailingOnly=TRUE)
infile = args[1] #input file
outdir = args[2] #output directory
nclusts = as.numeric(args[3]) #max number of clusters to consider
sim = F
if (length(args)>3){
   if (args[4] == 'sim'){ #simulated data
      sim = T}
 }


if (sim){ #we read simulated data from a tsv file
dat = read.table(infile, sep = '\t', header=T)
}else{ dat = readRDS( infile ) } #we read assignemnt data from an RDS file

numOfClusterPool = 1:nclusts
numOfRepeat = 1

t = as.numeric(Sys.time()) #measure time of calculation
results <- RunCcubePipeline(ssm = dat, 
                            numOfClusterPool = numOfClusterPool, 
                            numOfRepeat = numOfRepeat,
                            runAnalysis = T, 
                            runQC = T)
t = as.numeric(Sys.time())-t #store runtime
print(paste('runtime:', t))

summary = table(results$ssm$ccube_ccf_mean) #CCF and #mutations
cfs = as.numeric(names(summary))*100 #extract CCF data
sizes = as.numeric(summary) #extract #mutations

use_python('/local/data/public/kird/anaconda3/bin/python', required=T)
source_python('plot_data.py') #import plotting function
plot_data(cfs, sizes, dirname = outdir, sim=sim) #plot data

