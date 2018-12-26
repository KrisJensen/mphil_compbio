#questions
#can we assume chrom sizes in directory? yes; might want to highlight name of file
#rnw vs rmd - can use rnw (specify at top how to compile)
#how to sample? we decide - do something sensible and justify
#scale and see how close we get to integers
#do we not need noise in the un-copied segments?

#chr_sizes = read.table('chrom_sizes.txt')
#ends = cumsum(as.numeric(chr_sizes[,2]))



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


set.seed(26110940)

cn_change = function(results, tot.size, n_change, noise){
  #segs represents boundary segments (i.e. chromosomes or copy number regions)

  ps = cumsum( c(0, results[,'length']*results[,'cn_true']) / tot.size )
  pos = runif(1)
  seg = sum (as.numeric(pos > ps)); chr = as.character(results[seg, 'chr']);
  L = results[seg, 'length']
  
  site = round(runif(1)*L,0); size = round(runif(1, 10^6, 10^8), 0)
  init = max(0, site-floor(size/2)) #cannot cross segment boundary
  fin = min(L, site+ceiling(size/2)) #cannot cross segment boundary
  
  change = sample( c( -n_change:-1 , 1:n_change ), 1)
  new_cn = results[seg,'cn_true']+change #old copy number plus change
  new_cn = max(0, new_cn) # don't allow negative copy numbers
  
  if (noise){new_cns = c(new_cn, new_cn+rnorm(1,sd=0.10)) #add noise to observed
  }else{new_cns = c(new_cn, new_cn)} #obs and true are identical
  
  tot.size = tot.size + change*(fin-init) #update our total amount of DNA
  
  #cat('\nchr: ', chr, 'site: ',site, 'init: ', init, 'fin: ', fin,
  #    'size: ', size, 'real size: ', fin-init, '\n')
  #print(results[seg,])
  
  if ( isTRUE(all.equal( c(init,fin), c(0,L) ) )){ #complete overlap, update copy number
    results[seg,c('cn_true','cn_obs')]=new_cns
  }else if (init == 0){ #creates new CN at start of segment
    results[seg,'length'] = L - fin
    if(seg==1){results = rbind(c(NA, fin, new_cns), results)
    }else{results = rbind(results[1:(seg-1),], c(NA, fin, new_cns), results[-(1:(seg-1)),])}
    results[seg, 'chr'] = chr
  }else if (fin == L){ #new CN at end of segment
    results[seg,'length'] = init
    results = rbind(results[1:seg,], c(NA, L-init, new_cns), results[-(1:seg),])
    results[seg+1, 'chr'] = chr
  }else{ #falls within segment
    if(seg == 1){results = rbind(results[1,], c(NA, size, new_cns), results)
    }else{results = rbind(results[1:seg,], c(NA, size, new_cns), results[-(1:(seg-1)),])}
    results[seg,'length']= init
    results[seg+2,'length']= L-fin
    results[seg+1, 'chr'] = chr
  }
  #print(results)
  return(list(results, tot.size))
}

cnv_sim = function(N, n_change = 1, noise = 0, ploidy=2){
  #N is number of sequential events
  results = read.table('chrom_sizes.txt', stringsAsFactors = 0)
  results[,3:4] = rep(ploidy,24)
  colnames(results) = c('chr', 'length', 'cn_true', 'cn_obs') 
  tot.size = ploidy*3095677412
  for (i in 1:N){
    #print(results)
    new_res = cn_change(results, tot.size, n_change, noise)
    results = new_res[[1]]; tot.size = new_res[[2]]
  }
  #print(results)
  #return observed values
  #results = results[, c(1,2,4)]; colnames(results)=c('chr', 'length', 'cn') 
  results = results[, c(1,2,3)]; colnames(results)=c('chr', 'length', 'cn')
  if (noise){ results[,3] = results[,3] + rnorm(dim(results)[1], sd=0.1) }
  return(results)
}


plot_results = function(results, ploidy='NA', n_change='NA', noise=2, main=''){
  cum = cumsum(results[,2])
  ymax = min( max(results[,3], na.rm=T)+1, median(results[,3], na.rm=1)+5)
  plot(c(1, 3095677412), c(2,2), type='n', lty=2,
       xlim = c(1, 3095677412), ylim = c(-0.5, ymax),
       xlab = '', ylab = '', main = main, xaxt = 'n'
       )
  #print(cum)
  #print(results[,3])
  segments( c(0, cum[-length(cum)]), results[, 3],
            x1 = cum, col='blue', y1 = results[, 3], lwd=3  )
  chr_sizes = read.table('chrom_sizes.txt')
  cum = cumsum(as.numeric(chr_sizes[,2])); cum = cum[-length(cum)]
  segments( cum, -0.5, x1 = cum, ymax, lwd=0.5)
}

sim_more_cnvs = function(cnvs = c(3,10,25,100,200), n_change=1, noise=0, ploidy=2){
  results = vector('list', 0)
  for (cnv in cnvs){
    results[[cnv]] = cnv_sim(cnv, n_change, noise, ploidy)
    plot_results(results[[cnv]], ploidy, n_change, noise)
  }
  return(results)
}

check_results = function(results_df){
  #results[,4] = results[,3]; colnames(results) = c('chr', 'start', 'end', 'cn')
  chr_sizes = read.table('chrom_sizes.txt')
  for (i in 1:24){#check that sizes of chunks add up
    cat('\n', as.character(chr_sizes[i,1]), 'real', as.character(chr_sizes[i,2]),
        'sum', as.character(sum(results_df[results_df[,'chr']==chr_sizes[i,1],2])), '\n')
    #coords = cumsum(results[results[,'chr']==chr_sizes[i,1],2])
    #results[results[,'chr']==chr_sizes[i,1],2:3] =
    #              cbind( c(0, coords[-length(coords)]), coords )
    #print(results)
  }
}

clonal_proximity = function(results, n=2){
  L = sum(results[,'length'])
  #normalized metric runs from 0 if all is integer to 1 if all is at .5
  
  #d = sum(abs((results[,'cn'] - round(results[,'cn'],0)))*results[,'length'])  *  2/L
  d = ( sum( (results[,'cn'] - round(results[,'cn'],0) )^n * results[,'length'] )*4/L )^(1/n)
  
  sim = 1-d;# print(sim);
  return(sim)
}

normalize = function(results){
  results[,'cn']=results[,'cn']/median(results[,'cn'])
  return(results)}

test_proximity = function( results, scale_factors = 'default', main='default' ){
  res = results[ !is.na(results[,'cn']), ]
  sims = c()
  #print(min(res[,'cn'])); print(max(res[,'cn']))
  if (scale_factors == 'default'){
    #shifts = seq( -min(res[,'cn'])-1, #start screena at CN of -1 for lowest CN
    #              max(res[,'cn'])-min(res[,'cn']), #screen full range
    #              l = 100)
    scale_factors = seq( 0.5, 5, l=1000 ) #go from half-ploid to pentaploid reference
  }
  #print(scale_factors[1]); print(scale_factors[20])
  test_res=res
  for (scale in scale_factors){
    test_res[,'cn'] = res[,'cn']*scale
    sims = c(sims, clonal_proximity(test_res))
    #plot_results(test_res); print(shift, sims[length(sims)])
  }
  
  par(mfrow = c(3,1), mar = c(0.5,3.5,1.5,1))
  plot_results(results, main = paste('Relative Copy Number Profile', main) )
  
  simmax = max(sims); scalemax = scale_factors[which(sims==simmax)]
  results[,'cn'] = results[,'cn']*scalemax[1]
  
  plot_results(results, main = paste('Inferred Clonal Profile', main) )
  
  par(mar = c(3.5,3.5,1.5,1))
  plot(scale_factors, sims, type='l', main=paste('Scaling Profile', main),
       xlab = '', ylab = '', ylim = c(0,1))
  title(ylab="clonal similarity", line=2, cex.lab=1.2)
  title(xlab="scaling factor", line=2, cex.lab=1.2)
  
  
  return(list('results'=results, 'scaling'=scalemax, 'clonal'=simmax))
}

par(mfrow = c(5,1), mar = c(0.5,2.5,0.5,1), xaxs = 'i', yaxs = 'i', lend=1)
res1 = sim_more_cnvs( noise=0, n_change=1)
res2 = sim_more_cnvs( noise=1, n_change=2)
res3 = sim_more_cnvs( noise=1, n_change=2, ploidy=4)

res1_norm = normalize(res1[[10]])
shift1 = test_proximity(res1_norm, main = 'Without Noise')


res4 = cnv_sim (10, n_change = 1, noise = 1, ploidy=2)
res4_norm = normalize(res4)
shift4 = test_proximity(res4_norm, main = 'With Gaussian Noise')


parse_segs = function(segs){
  newsegs = list()
  chr_sizes = read.table('chrom_sizes.txt')
  for (i in 1:length(segs)){
    seg = segs[[i]]
    seg[,'start'] = as.numeric(seg[,'start'])
    seg[,'end'] = as.numeric(seg[,'end'])
    seg[,'cn'] = as.numeric(seg[,'cn'])
    for (j in c(1:22,'X','Y')){
      seg[seg[,'chr'] == j, 'chr'] = paste0('chr',j)
    }
    newseg = data.frame('chr' = c(seg[1,'chr']),
                        'length' = c(seg[1,'end'] - seg[1,'start']+1),
                        'cn' = c(seg[1,'cn']), stringsAsFactors = FALSE)
    #print(newseg)
    if (seg[1,'start'] != 1){
      newseg = rbind(list(seg[1,'chr'], seg[1,'start']-1,  NA), newseg)
    }
    if ((seg[1,'end']+1 != seg[2,'start']) & (seg[1, 'chr'] == seg[2, 'chr'])){
      newseg = rbind(newseg, list(seg[1,'chr'],
                                  seg[2,'start'] - seg[1,'end']-1,  NA))
    }
    #print(newseg)
    for (n in 2:(dim(seg)[1]-1)){
      #print(seg[n,'chr'])
      if (seg[n, 'chr'] != seg[n-1, 'chr']){
        if (seg[n,'start'] != 1){
          newseg = rbind(newseg, list(seg[n,'chr'], seg[n,'start']-1,  NA))
        }
      }
      newseg = rbind(newseg, list(seg[n,'chr'], seg[n,'end'] - seg[n,'start']+1,  seg[n,'cn']))
      if ((seg[n,'end']+1 != seg[n+1,'start']) & (seg[n, 'chr'] == seg[n+1, 'chr'])){
        newseg = rbind(newseg, list(seg[n,'chr'],
                        seg[n+1,'start'] - seg[n,'end']-1,  NA))
      }
      if (seg[n, 'chr'] != seg[n+1, 'chr']){
        if (seg[n,'end'] != chr_sizes[chr_sizes[,1] == seg[n,'chr'],2]){
          newseg = rbind(newseg, list(seg[n,'chr'], 
                                  chr_sizes[chr_sizes[,1] == seg[n,'chr'],2]-seg[n,'end'],NA))
        }
      }
    }
    newseg = rbind(newseg, list(seg[dim(seg)[1],'chr'],
                                seg[dim(seg)[1],'end'] - seg[dim(seg)[1],'start']+1,
                                seg[dim(seg)[1],'cn']))
    if (seg[dim(seg)[1],'end'] != chr_sizes[chr_sizes[,1] == seg[n,'chr'],2]){
      newseg = rbind(newseg, list(seg[dim(seg)[1],'chr'], 
                                  chr_sizes[chr_sizes[,1] == seg[dim(seg)[1],'chr'],2]
                                  -seg[dim(seg)[1],'end'], NA))
    }
    newsegs[[i]] = newseg
  }
  names(newsegs) = names(segs)
  return(newsegs)
  #put NAs where no info
}
length_to_seg = function(results){
  segs = data.frame('chr' = c(results[1,'chr']), 'start' = c(1),
                    'end' = c(results[1,'length']), 'cn' = c(results[1,'cn']),
                    stringsAsFactors = 0)
  x = results[1,'length']
  
  for (i in 2:dim(results)[1]){
    if (results[i, 'chr'] == results[i-1, 'chr']){
      segs[i,] = list(results[i,'chr'], x+1, x+results[i,'length'], results[i,'cn'])
      x = x + results[i, 'length']
    }else{
      segs[i,] = list(results[i,'chr'], 1, results[i,'length'], results[i,'cn'])
      x = results[i, 'length']
    }
  }
  segs[,'cn'] = round(segs[,'cn'],2)
  return(segs)
}

# segs = readRDS('relative_segment_tables.rds')
# news = parse_segs(segs)
# plot_results(news[[1]])
#
# for (i in 1:5){
#   plot_results(news[[i]], main = paste0('Absolute profile ', names(news)[i]) )
#   scale_rds = test_proximity(news[[i]], main = paste0('Scaling profile ', names(news)[i]));
#   print(scale_rds[[2]])
#   plot_results(scale_rds[[1]], main = paste0('Absolute profile ', names(news)[i]))
#   back = length_to_seg(scale_rds[[1]])
# }
