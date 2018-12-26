#read in mmseq output files and construct DESeq dataset

library(DESeq2)

#use mmseq2counts.R courtesy of Ernest Turro (RNA seq practical)
source("mmseq2counts.R")

#sample names/output files
outs = c('adipose_53', 'spleen_51', 'pancreas_51', 'colon_51', 'heart_51', 'breast-epithelium_51',
	'liver_emb', 'liver_6', 'liver_53',
	'cardiac_1', 'cardiac_2', 'epithelium_1', 'epithelium_2',
	'hepatocyte_1', 'hepatocyte_2', 'neural-progenitor_1',
	'neural-progenitor_2')#, 'Mmus_1', 'Mmus_2')

fileouts = paste0('../bam/cdna/',outs,'_mm/',outs,'.mmseq')

#parse mmseq utputfiles
tab_iso = mmseq2counts(fileouts, sample_names = outs)

labs = c('adipose', 'spleen', 'pancreas', 'colon', 'heart', 'breast-epithelium',
        'liver_embryo', 'liver_6', 'liver_53',
        'cardiac', 'cardiac', 'epithelium', 'epithelium',
        'hepatocyte', 'hepatocyte', 'neural-progenitor',
        'neural-progenitor')#, 'Mmus', 'Mmus')

#store results
saveRDS(tab_iso, file = "tab_isoforms.rds")

#construct DESeq dataset
coldata = matrix( labs, length(outs), 1); rownames(coldata) = outs; colnames(coldata) = 'Tissue'
cds_iso = DESeqDataSetFromMatrix(round(tab_iso), colData = coldata, design = ~ Tissue )
cds_iso = estimateSizeFactors( cds_iso )

#save DESeq dataset
saveRDS(cds_iso, file = "dds_isoforms.rds")

