#parse mmseq data for liver samples from different species

library(DESeq2)

#courtesy of Ernesto Turro from RNAseq practical
source("mmseq2counts.R")

outs = c('liver_emb', 'liver_6', 'liver_53',
	'hepatocyte_1', 'hepatocyte_2', 'Mmus_1', 'Mmus_2')

fileouts = paste0('../bam/cdna/',outs,'_mm/',outs,'.mmseq')

#parse mmseq output files
tab_iso = mmseq2counts(fileouts, sample_names = outs)

labs = c('liver_embryo', 'liver_6', 'liver_53',
        'hepatocyte', 'hepatocyte', 'Mmus', 'Mmus')

saveRDS(tab_iso, file = "tab_spec_isoforms.rds")

#analyze human and mouse data in parallel
coldata1 = matrix( labs[1:5], 5, 1); rownames(coldata1) = outs[1:5]; colnames(coldata1) = 'Tissue'
coldata2 = matrix( labs[6:7], 2, 1); rownames(coldata2) = outs[6:7]; colnames(coldata2) = 'Tissue'

cds_iso1 = DESeqDataSetFromMatrix(round(tab_iso)[,1:5], colData = coldata1, design = ~ Tissue )
cds_iso1 = estimateSizeFactors( cds_iso1 )

cds_iso2 = DESeqDataSetFromMatrix(round(tab_iso)[,6:7], colData = coldata2, design = ~ 1 )
cds_iso2 = estimateSizeFactors( cds_iso2 )

#combine human and mouse data
counts = cbind(counts(cds_iso1, normalized=T), counts(cds_iso2, normalized=T))

#save data as count matrix
saveRDS(tab_iso, file = "tab_spec_isoforms.rds")

#make DEseq dataset
cds_iso = DESeqDataSetFromMatrix(round(counts), colData = rbind(coldata1, coldata2), design = ~ Tissue )
saveRDS(cds_iso, file = "dds_spec_isoforms.rds")

