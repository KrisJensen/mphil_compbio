require('ccube')
require('dplyr')

numSnv <- 500
ccfSet <- c(1, 0.4, 0.6) # true ccf pool
ccfTrue <- sample(ccfSet, numSnv, c(0.5,0.2,0.3), replace = T) # simulate true clusters
purity <- 0.9
cnPoolMaj <- c(1,2,3,4) # a pool of possible major copy numbers
cnPoolMin <- c(0,1,2) # a pool of possible minor copy numbers
cnPoolMajFractions <- c(0.30, 0.30, 0.2,0.2) # prevalence of possible major copy numbers
cnPoolMinFractions <- c(1/4, 1/2, 1/4) # prevalence of possible minor copy numbers
cnProfile = GenerateCopyNumberProfile(cnPoolMaj, cnPoolMin, 
                                      cnPoolMajFractions, cnPoolMinFractions, numSnv)
head(cnProfile) # column 1: minor copy number, column 2: major copy number, column 3: total copy number 


baseDepth = 50
mydata <- data.frame(mutation_id = paste0("ss","_", seq_len(numSnv)) ,
                     ccf_true = ccfTrue,
                     minor_cn = cnProfile[,1],
                     major_cn = cnProfile[,2],
                     total_cn = cnProfile[,3], 
                     purity = purity,
                     normal_cn = 2)
mydata <- dplyr::mutate(rowwise(mydata),
                        mult_true = sample(seq(1,if (major_cn ==1) { 1 } else {major_cn}), 1), # simulate multiplicity
                        vaf = cp2ap(ccf_true, purity, normal_cn, total_cn, total_cn, mult_true), # simulate vaf
                        total_counts = rpois(1, total_cn/2 * baseDepth), # simulate total read counts
                        var_counts = rbinom(1, total_counts, vaf),  # simulate variant read counts
                        ref_counts = total_counts - var_counts)
head(mydata)


