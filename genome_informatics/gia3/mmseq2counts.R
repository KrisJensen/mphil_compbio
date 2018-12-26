mmseq2counts <- function( mmseq_files = grep("gene|identical", dir( pattern="\\.mmseq$" ), value=TRUE, invert=TRUE),
		 sample_names = sub("\\.mmseq","", mmseq_files) ) {
	
	cat( "Using\n\tmmseq files:", mmseq_files, "\n\tsample names:", sample_names, "\n" )
	if( length(mmseq_files) == 0 || length(sample_names)  != length(mmseq_files) ) {
		cat( "ERROR:\n\tplease define mmseq files and sample names for each file.\n" )
		return(NULL)
	}
	
	# read in mmseq data
	data = list()
	lengths = list()
  n_reads <- c()
	cat( "Reading mmseq output files\n" )
	for( i in 1:length(mmseq_files) ) {
		cat( "\t", mmseq_files[i], "\t" )
		dat = read.table( mmseq_files[i], header=TRUE, stringsAsFactors=FALSE )
    if(any(colnames(dat)=="log_mu")) {
      data[[sample_names[i]]] = dat$log_mu
    } else {
      data[[sample_names[i]]] = dat$mean_log_mu_gibbs
    }
	print(colnames(dat))
		names(data[[sample_names[i]]]) = dat[,grepl("feature_id", colnames(dat))]
    lengths[[sample_names[i]]] = dat[,grepl("true_length", colnames(dat))]
		names(lengths[[sample_names[i]]]) = names(data[[sample_names[i]]])
    ta <- readLines(mmseq_files[i], n=100) # header should be no more than 100 lines
    ta <- strsplit(ta[grep("# Mapped fragments", ta)], " ")[[1]]
    n_reads <- c(n_reads, as.numeric(ta[length(ta)]))
		cat( "found", n_reads[length(n_reads)], "aligned reads (or read pairs)\n" )
	}
  names(n_reads) <- sample_names
	
	cat( "Sorting transcripts\n" )
	allnames = sort(unique(unlist(lapply( data, function(x) names(x) ))))
	tab = data.frame( row.names=allnames )
	for( e in data ) {
		tab = cbind(tab, rep(NA, nrow(tab)))
		tab[match(names(e), allnames),ncol(tab)] = e
	}
	colnames(tab) = names(data)
	t_lengths = data.frame( row.names=allnames )
	for( e in lengths ) {
		t_lengths = cbind(t_lengths, rep(NA, nrow(t_lengths)))
		t_lengths[match(names(e), allnames),ncol(t_lengths)] = e
	}
	colnames(t_lengths) = names(data)
	
	# un-log the estimates 
	cat( "Unlogging and unstandardising the estimates\n" )
	tab = exp(1)^(tab)
	
	# multiply by library size
	utab = sapply( colnames(tab), function(x) {	tab[,x]*n_reads[x] })
	rownames(utab)=rownames(tab)
	
	# multiply by transcript length and divide by 10^9 - i.e. unnormalize
	utab = sapply( colnames(tab), function(x) { utab[,x] = (utab[,x]*t_lengths[,x])/(1e3*1e6) } )
	utab[is.na(utab)] = 0 # the ones that had no expression in MMSEQ
	
	return(utab)
}
