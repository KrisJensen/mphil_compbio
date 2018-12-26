#parse bam files aligned to genome for Homo sapiens

library(biomaRt)
library(GenomicAlignments)
library(DESeq2)

ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
annot<-getBM(
	c("chromosome_name",
	"strand",
	"ensembl_gene_id",
	"ensembl_transcript_id",
	"ensembl_exon_id",
	"start_position", "end_position",
	"transcript_start",
	"transcript_end",
	"transcript_length",
	"exon_chrom_start",
	"exon_chrom_end"),
	mart=ensembl,
	filters="chromosome_name",
	values=c("13")
	)


#samples of interest
sams = c('adipose_53', 'spleen_51', 'pancreas_51', 'colon_51', 'heart_51', 'breast-epithelium_51',
        'liver_emb', 'liver_6', 'liver_53',
        'cardiac_1', 'cardiac_2', 'epithelium_1', 'epithelium_2',
        'hepatocyte_1', 'hepatocyte_2', 'neural-progenitor_1',
        'neural-progenitor_2')#, 'Mmus_1', 'Mmus_2')

#labels
labs = c('adipose', 'spleen', 'pancreas', 'colon', 'heart', 'breast-epithelium',
        'liver_embryo', 'liver_6', 'liver_53',
        'cardiac', 'cardiac', 'epithelium', 'epithelium',
        'hepatocyte', 'hepatocyte', 'neural-progenitor',
        'neural-progenitor')#, 'Mmus', 'Mmus')

geneCounts = c()
for(s in sams) {
	cat("Reading in ", s, "\n")
	#read bam file
	alnRanges <- readGAlignments(paste0("../bam/dna/", s, "_th/accepted_hits.bam"))
	exonRanges <- GRanges(seqnames = Rle(annot$chromosome_name),
	ranges = IRanges( start=annot$exon_chrom_start,
		end=annot$exon_chrom_end),
		strand = Rle( annot$strand ),
		exon=annot$ensembl_exon_id,
		gene=annot$ensembl_gene_id )
	strand(exonRanges) <- "*"
	exonCounts <- countOverlaps( exonRanges, alnRanges ) #assign to gene
	names( exonCounts ) <- elementMetadata(exonRanges)$gene
	splitCounts <- split(exonCounts, names(exonCounts) )
	geneCounts <- cbind(geneCounts, sapply( splitCounts, function(x) sum(x) ) )
}
colnames(geneCounts) <- sams #each column is a sample
head(geneCounts)

saveRDS(geneCounts, file = "matrix_genes.rds") #store for differential expression analysis

coldata = matrix( labs, length(sams), 1)
rownames(coldata) = colnames(geneCounts)
colnames(coldata) = 'Tissue'
coldata = data.frame(coldata)
#construct read dataset
dds <- DESeqDataSetFromMatrix(geneCounts, colData = coldata, design = ~ Tissue)
dds <- estimateSizeFactors(dds)
cat(sizeFactors(dds))

saveRDS(dds, file = "dds_genes.rds") #store for plotting

