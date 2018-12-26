#questions
#can we assume chrom sizes in directory? yes; might want to highlight name of file
#rnw vs rmd - can use rnw (specify at top how to compile)
#how to sample? we decide - do something sensible and justify
#scale and see how close we get to integers
#do we not need noise in the un-copied segments?

# Considerations: CN changes happen uniformly along the genome; spread out symmetrically
# from point of event (equivalent to randomly drawing a second break point at either side)
# Proportional to the copy number of the region
# this also means that once a region is completely lost, this loss is irreversible
# and we do not observe CN changes \textit{from} zero.
# When introducing noise, we assume that there is still an underlying discrete
# ground truth biology whereas noise occurs from errors in measurements or changes in
# transcriptional activity from translocations or other such effects.
# We thus let the ground truth biology dictate probabilities of CN changes while
# adding noise to our observed output. (cn_obs, cn_true in df)


set.seed(26110940) #set seed for repeatable simulations

cn_change = function(results, tot.size, n_change, noise){
  #simulates a single copy number change given a data frame of segments (results)
  #a cn change is initiated uniformly at random across the genome
  #it then extends symmetrically in both directions with a length
  #between 1 and 100 Mb picked uniformly at random, bounded by the segment in which it occurs

  #probability of initating a cn change at each basepair
  #is proportional to its current copy number
  ps = cumsum( c(0, results[,'length']*results[,'cn']) / tot.size )
  pos = runif(1)
  #first find the segment in which a CN change occurs, then the actual point of change
  seg = sum (as.numeric(pos > ps)); chr = as.character(results[seg, 'chr']);
  L = results[seg, 'length']
  site = round(runif(1)*L,0); size = round(runif(1, 10^6, 10^8), 0)
  
  init = max(0, site-floor(size/2)) #cn change cannot go past segment start
  fin = min(L, site+ceiling(size/2)) #cn change cannot go past segment end
  
  change = sample( c( -n_change:-1 , 1:n_change ), 1)
  new_cn = results[seg,'cn']+change #old copy number plus change
  new_cn = max(0, new_cn) # don't allow negative copy numbers
  
  tot.size = tot.size + change*(fin-init) #update our total amount of DNA
  
  if ( isTRUE(all.equal( c(init,fin), c(0,L) ) )){ #complete overlap, update copy number
    results[seg,'cn']=new_cn
  }else if (init == 0){ #creates new CN at start of segment
    results[seg,'length'] = L - fin
    if(seg==1){results = rbind(list(chr, fin, new_cn), results)
    }else{results = rbind(results[1:(seg-1),], list(chr, fin, new_cn), results[-(1:(seg-1)),])}
  }else if (fin == L){ #new CN at end of segment
    results[seg,'length'] = init
    results = rbind(results[1:seg,], list(chr, L-init, new_cn), results[-(1:seg),])
  }else{ #cn change falls within segment
    if(seg == 1){results = rbind(results[1,], list(chr, size, new_cn), results)
    }else{results = rbind(results[1:seg,], list(chr, size, new_cn), results[-(1:(seg-1)),])}
    results[seg,'length']= init
    results[seg+2,'length']= L-fin
  }

  return(list(results, tot.size))
}

cnv_sim = function(N, n_change = 1, noise = 0, ploidy=2){
  #simulates N copy number changes in a hypothetical human genome
  #n_change is max copy number change in a single event
  #noise specifies whether we add Gaussian noise, ploidy specifies initial ploidy
  
  results = read.table('chrom_sizes.txt', stringsAsFactors = 0)
  results[,3] = rep(ploidy,24) #initiate genome
  colnames(results) = c('chr', 'length', 'cn') 
  tot.size = ploidy*3095677412 #total number of base pairs
  
  for (i in 1:N){
    new_res = cn_change(results, tot.size, n_change, noise) #simulate cn change
    results = new_res[[1]]; tot.size = new_res[[2]]
  }

  colnames(results)=c('chr', 'length', 'cn')
  #Add Gaussian noise to all segments if specified
  if (noise){ results[,3] = results[,3] + rnorm(dim(results)[1], sd=0.1) }
  return(results)
}


plot_results = function(results, ploidy='NA', n_change='NA', noise=2, main=''){
  #given a data frame of segment lengths and copy numbers, plots a copy number profile
  cum = cumsum(results[,2]) #x boundaries
  ymax = min( max(results[,3], na.rm=T)+1, median(results[,3], na.rm=1)+5)
  plot(c(1, 3095677412), c(2,2), type='n', lty=2,
       xlim = c(1, 3095677412), ylim = c(-0.5, ymax),
       xlab = '', ylab = '', main = main, xaxt = 'n'
       ) #initiatie plot

  segments( c(0, cum[-length(cum)]), results[, 3], #plot copy numbers
            x1 = cum, col='blue', y1 = results[, 3], lwd=3  )
  chr_sizes = read.table('chrom_sizes.txt')
  cum = cumsum(as.numeric(chr_sizes[,2])); cum = cum[-length(cum)]
  segments( cum, -0.5, x1 = cum, ymax, lwd=0.5) #plot chromosome delineations
}

sim_more_cnvs = function(cnvs = c(3,10,25,100,200), n_change=1, noise=0, ploidy=2){
  #function for simulating multiple CNVs and putting in a plot
  par(mfrow = c(length(cnvs),1), mar = c(0.5,2.5,0.5,1), xaxs = 'i', yaxs = 'i', lend=1)
  results = vector('list', 0)
  for (cnv in cnvs){
    results[[cnv]] = cnv_sim(cnv, n_change, noise, ploidy) #simulate changes
    plot_results(results[[cnv]], ploidy, n_change, noise) #plot result
  }
  return(results)
}

check_results = function(results_df){
  #given a set data frame of segments, checks that
  #the length of all segments in a chromosome adds up to the total chromsome length
  chr_sizes = read.table('chrom_sizes.txt')
  for (i in 1:24){#check that sizes of chunks add up
    cat('\n', as.character(chr_sizes[i,1]), 'real', as.character(chr_sizes[i,2]),
        'sum', as.character(sum(results_df[results_df[,'chr']==chr_sizes[i,1],2])), '\n')
  }
}

clonal_proximity = function(results, n=2){
  #quantifies how close a copy number profile is to being clonal
  #we use a root mean square distance metric normalized by segment length
  #this is 0 if all segments have integer cn, 1 if they are all half-integer
  
  L = sum(results[,'length'])
  #d = sum(abs((results[,'cn'] - round(results[,'cn'],0)))*results[,'length'])  *  2/L
  d = ( sum( (results[,'cn'] - round(results[,'cn'],0) )^n * results[,'length'] )*4/L )^(1/n)
  sim = 1-d; #similarity is one minus normalized distance
  return(sim)
}

normalize = function(results){
  #performs median normalization
  results[,'cn']=results[,'cn']/median(results[,'cn'])
  return(results)}

test_proximity = function( results, scale_factors = 'default', main='default' ){
  #given a dataframe of relative copy numbers (results)
  #finds a scaling to convert it to absolute copy number using the
  #clonal_proximity() function. Plots relative and absolute copy number profiles
  #together with the clonal similarity as a function of scaling
  res = results[ !is.na(results[,'cn']), ]
  sims = c()

  if (scale_factors == 'default'){
    scale_factors = seq( 0.5, 5, l=1000 ) #go from haploid to pentaploid reference
  }

  test_res=res
  for (scale in scale_factors){
    test_res[,'cn'] = res[,'cn']*scale #scale copy numbers
    sims = c(sims, clonal_proximity(test_res)) #assess the quality of our result
  }
  
  par(mfrow = c(3,1), mar = c(0.5,3.5,1.5,1)) #plot relative cn profile
  plot_results(results, main = paste('Relative Copy Number Profile', main) )
  
  simmax = max(sims); scalemax = scale_factors[which(sims==simmax)]
  results[,'cn'] = results[,'cn']*scalemax[1] #calculate absolute profile
  
  #plot aboslute profile
  plot_results(results, main = paste('Inferred Clonal Profile', main) )
  
  par(mar = c(3.5,3.5,1.5,1)) #plot clonal similarity vs scaling
  plot(scale_factors, sims, type='l', main=paste('Scaling Profile', main),
       xlab = '', ylab = '', ylim = c(0,1))
  title(ylab="clonal similarity", line=2, cex.lab=1.2)
  title(xlab="scaling factor", line=2, cex.lab=1.2)
  
  return(list('results'=results, 'scaling'=scalemax, 'clonal'=simmax))
}

# res1 = sim_more_cnvs( noise=0, n_change=1)
# res2 = sim_more_cnvs( noise=1, n_change=2)
# res3 = sim_more_cnvs( noise=1, n_change=2, ploidy=4)
# 
# res1_norm = normalize(res1[[10]])
# shift1 = test_proximity(res1_norm, main = 'Without Noise')
# 
# res3_norm = normalize(res3[[10]])
# shift3 = test_proximity(res3_norm, main = 'Tetraploid')
# 
# res4 = cnv_sim(10, n_change = 1, noise = 1, ploidy=2)
# res4_norm = normalize(res4)
# shift4 = test_proximity(res4_norm, main = 'With Gaussian Noise')
# 
# 
parse_segs = function(segs){
  #given a list of dataframes with structure [chr, start, end, cn]
  #returns a list of dataframes with structure [chr, length, cn]
  #fills in segments for which we don't have data with NA

  newsegs = list()
  chr_sizes = read.table('chrom_sizes.txt')
  for (i in 1:length(segs)){
    seg = segs[[i]]
    seg[,'start'] = as.numeric(seg[,'start'])
    seg[,'end'] = as.numeric(seg[,'end'])
    seg[,'cn'] = as.numeric(seg[,'cn'])

    for (j in c(1:22,'X','Y')){ #standardize chromosome labels
      seg[seg[,'chr'] == j, 'chr'] = paste0('chr',j)
    }
    newseg = data.frame('chr' = c(seg[1,'chr']),
                        'length' = c(seg[1,'end'] - seg[1,'start']+1),
                        'cn' = c(seg[1,'cn']), stringsAsFactors = FALSE)
    if (seg[1,'start'] != 1){ #if we don't have data for beginning of chr1, add NA
      newseg = rbind(list(seg[1,'chr'], seg[1,'start']-1,  NA), newseg)
    }
    if ((seg[1,'end']+1 != seg[2,'start']) & (seg[1, 'chr'] == seg[2, 'chr'])){
      #if we have a gap between segments, add NA
      newseg = rbind(newseg, list(seg[1,'chr'],
                                  seg[2,'start'] - seg[1,'end']-1,  NA))
    }
    for (n in 2:(dim(seg)[1]-1)){ #run through segments
      if (seg[n, 'chr'] != seg[n-1, 'chr']){
        #if new chromsome
        if (seg[n,'start'] != 1){ #if we don't have data for beginning of chr, add NA
          newseg = rbind(newseg, list(seg[n,'chr'], seg[n,'start']-1,  NA))
        }
      }
      newseg = rbind(newseg, list(seg[n,'chr'], seg[n,'end'] - seg[n,'start']+1,  seg[n,'cn']))
      if ((seg[n,'end']+1 != seg[n+1,'start']) & (seg[n, 'chr'] == seg[n+1, 'chr'])){
        #if we have a gap between segments, add NA
        newseg = rbind(newseg, list(seg[n,'chr'],
                        seg[n+1,'start'] - seg[n,'end']-1,  NA))
      }
      if (seg[n, 'chr'] != seg[n+1, 'chr']){
        #if we don't have data for end of chr, add NA
        if (seg[n,'end'] != chr_sizes[chr_sizes[,1] == seg[n,'chr'],2]){
          newseg = rbind(newseg, list(seg[n,'chr'],
                                  chr_sizes[chr_sizes[,1] == seg[n,'chr'],2]-seg[n,'end'],NA))
        }
      }
    }
    #fill out the rest of the final chromosome
    newseg = rbind(newseg, list(seg[dim(seg)[1],'chr'],
                                seg[dim(seg)[1],'end'] - seg[dim(seg)[1],'start']+1,
                                seg[dim(seg)[1],'cn']))
    if (seg[dim(seg)[1],'end'] != chr_sizes[chr_sizes[,1] == seg[n,'chr'],2]){
      newseg = rbind(newseg, list(seg[dim(seg)[1],'chr'],
                                  chr_sizes[chr_sizes[,1] == seg[dim(seg)[1],'chr'],2]
                                  -seg[dim(seg)[1],'end'], NA))
    }
    newsegs[[i]] = newseg #add to list of segments with new format
  }
  names(newsegs) = names(segs)
  return(newsegs)
}
length_to_seg = function(results){
  #given a dataframe with structure [chr, length, cn]
  #returns a dataframes with structure [chr, start, end, cn]
  segs = data.frame('chr' = c(results[1,'chr']), 'start' = c(1),
                    'end' = c(results[1,'length']), 'cn' = c(results[1,'cn']),
                    stringsAsFactors = 0)
  x = results[1,'length'] #use x as positional counter along chromsome

  for (i in 2:dim(results)[1]){
    if (results[i, 'chr'] == results[i-1, 'chr']){ #if same chromosome
      segs[i,] = list(results[i,'chr'], x+1, x+results[i,'length'], results[i,'cn'])
      x = x + results[i, 'length']
    }else{ #if starting new chromsome, reset counter
      segs[i,] = list(results[i,'chr'], 1, results[i,'length'], results[i,'cn'])
      x = results[i, 'length']
    }
  }
  return(segs)
}

 segs = readRDS('relative_segment_tables.rds')
 news = parse_segs(segs)
# 
# for (i in 1:5){
#   scale_rds = test_proximity(news[[i]], main = names(news)[i]);
#   print(scale_rds[2:3])
#   back = length_to_seg(scale_rds[[1]])
# }

test_proximity_ss = function( results, scale_factors = 'default',
                              shift_factors = 'default', main='default' ){
  res = results[ !is.na(results[,'cn']), ]
  
  #print(min(res[,'cn'])); print(max(res[,'cn']))
  if (scale_factors == 'default'){
    #shifts = seq( -min(res[,'cn'])-1, #start screena at CN of -1 for lowest CN
    #              max(res[,'cn'])-min(res[,'cn']), #screen full range
    #              l = 100)
    scale_factors = seq( 0.5, 5, l=91 ) #go from half-ploid to pentaploid reference
  }
  if (shift_factors == 'default'){
    shift_factors = seq( -0.3,0.05, l=101)
  }
  sims = matrix(,length(shift_factors), length(scale_factors))
  plotsims = matrix(,length(scale_factors), length(shift_factors))
  #print(scale_factors[1]); print(scale_factors[20])
  test_res=res
  for (shift in 1:length(shift_factors)){
    for (scale in 1:length(scale_factors)){
      test_res[,'cn'] = (res[,'cn']+shift_factors[shift])*scale_factors[scale]
      sim = clonal_proximity(test_res)
      sims[shift,scale] = sim
      plotsims[scale, length(shift_factors)+1-shift] = sim
      #plot_results(test_res); print(shift, sims[length(sims)])
    }
  }
  
  par(mfrow = c(2,1), mar = c(0.5,3.5,1.5,1))
  plot_results(results, main = paste('Relative Copy Number Profile', main) )
  
  simmax = max(sims); param_max = which(sims==simmax, arr.ind=1)
  shiftmax = shift_factors[param_max[1,1]]; scalemax = scale_factors[param_max[1,2]]
  print(simmax)
  print(param_max)
  print(shiftmax)
  print(scalemax)
  results[,'cn'] = (results[,'cn']+shiftmax)*scalemax
  
  plot_results(results, main = paste('Inferred Clonal Profile', main) )
  
  par(mfrow = c(1,1), mar = c(4,4,4,4))
  image(plotsims, axes = 0, xlab='', ylab='')
  axis(3, at = seq(0, 1, length = 10),
       labels=seq(scale_factors[1], scale_factors[length(scale_factors)], l=10)
       ,srt=45,tick=TRUE)
  axis(2, at = seq(0, 1, length = 5),
       labels=seq(shift_factors[length(shift_factors)], shift_factors[1], l=5)
       ,srt=45,tick=TRUE)
  title(ylab="shift", line=2, cex.lab=1.2)
  title(xlab="scaling factor", line=0, cex.lab=1.2)
  
  return(list('results'=results, 'scaling'=scalemax, 'shift'=shiftmax, 'clonal'=simmax))
}

test_proximity_ss(news[[3]])
test_proximity_ss(news[[5]])