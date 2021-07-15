gettingSP_ldscMR = function(Data,EXP,OUT,do_ldsc=FALSE){
  library(TwoSampleMR)
  library(GenomicSEM)
  
  file_output = paste0(EXP,"-",OUT,"_MR-SP.csv")
  
  nX = mean(Data$n_complete_samples.x)  #Get sample size for trait X
  nY = mean(Data$n_complete_samples.y)  #Get sample size for trait Y
  
  bX = Data$tstat.x/sqrt(nX)   #Get standardised beta for trait X
  bY = Data$tstat.y/sqrt(nY)   #Get standardised beta for trait Y

  X = select(Data, rsid, chr, alt, ref, beta.x, se.x, pval.x, tstat.x) %>% rename(A1 = alt, A2 = ref, unstdb = beta.x, sderr = se.x, pval = pval.x, tstat=tstat.x)
  Y = select(Data, rsid, chr, alt, ref, beta.y, se.y, pval.y, tstat.y) %>% rename(A1 = alt, A2 = ref, unstdb = beta.y, sderr = se.y, pval = pval.y, tstat=tstat.y)
  X$beta = bX
  Y$beta = bY
  X$se = 1/sqrt(nX)
  Y$se = 1/sqrt(nY)
  
  if(do_ldsc){
    # slice, every 10th SNP for faster computation
    X %>%
      slice(seq(1, nrow(X), by=10)) -> X_filtered
    nrow(X_filtered)
    Y %>%
      slice(seq(1, nrow(Y), by=10)) -> Y_filtered
    nrow(Y_filtered)
    
    write.table(X_filtered, file = paste0(EXP, "_GWAS.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names=TRUE)
    write.table(Y_filtered, file = paste0(OUT, "_GWAS.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names=TRUE)
    ## First you need to "munge" (pre-process, align ...) the data:
    # the GWAS.tsv files should be in the folder
    munge( paste0(EXP, "_GWAS.txt"), 
           "/project/data/genomicSEM/data/w_hm3.noMHC.snplist",
           trait.names=EXP,
           N=nX) 
    munge( paste0(OUT, "_GWAS.txt"), 
           "/project/data/genomicSEM/data/w_hm3.noMHC.snplist",
           trait.names=OUT,
           N=nY) 
    traits <- c(paste0(EXP, ".sumstats.gz"),paste0(OUT, ".sumstats.gz"))
    sample.prev <- c(NA,NA) # continuous traits
    population.prev <- c(NA,NA) # continuous traits
    ld<-"/project/data/genomicSEM/data/eur_w_ld_chr/"
    wld <- "/project/data/genomicSEM/data/eur_w_ld_chr/"
    trait.names<-c("GWAS1", "GWAS2")
    LDSCoutput <- ldsc(traits, 
                       sample.prev, 
                       population.prev, 
                       ld, 
                       wld, 
                       trait.names)
    save(LDSCoutput, file="Pfactor.RData")
    
    #i_X = LDSCoutput$I[1,1]
    #i_Y = LDSCoutput$I[2,2]
    i_XY = LDSCoutput$I[1,2]
    #h2_x_ldsc = LDSCoutput$S[1,1]
    #h2_y_ldsc = LDSCoutput$S[2,2]
  } else {i_XY = runif(1,-0.1,0.1)}
  
  
  ### Standard MR
  
  ## Get significant SNPs above certain Z-statistic corresponding to set p-value
  prune_X = function(zX,p_limit=1e-5){
    z_limit=abs(qnorm(0.5*p_limit))
    ind_keep=which(abs(zX)>z_limit)
    ind_keep=unique(ind_keep)
    ind_keep=list(ind_keep)
    return(ind_keep)
  }
  
  ## Taken from Jonathan Sulc to create bins that fit a maximum of 50k SNPs (max threshold for clumping)
  snp_bin  =  function( snp_ranks,
                        chunk_size = SNP_CHUNK_SIZE ){
    if (nrow( snp_ranks ) == 0) {
      return()
    }
    
    max_chr  =  snp_ranks$chr %>%
      table %>%
      cumsum %>%
      (function(x) x < chunk_size) %>%
      (function(x) names(x)[ max(which(x)) ] ) %>%
      as.numeric
    if (is.na( max_chr )) {
      max_chr = min( snp_ranks$chr )
    }
    
    bin = snp_ranks %>%
      filter( chr <= max_chr ) %>%
      list
    return( c( bin,
               snp_bin( snp_ranks[ snp_ranks$chr > max_chr, ],
                        chunk_size ) ) )
  }
  
  # Set values
  pval=2*pnorm(-abs(5.45))
  pval1=2*pnorm(-abs(4))
  reverse_t_threshold  =  qnorm( 5e-2 )
  
  ### Forward
  mr_dataX = cbind.data.frame(SNP = X$rsid, beta = X$beta, se = X$se, effect_allele = X$A1, other_allele = X$A2, chr=X$chr, Phenotype=EXP, tstat=X$tstat )
  mr_dataY = cbind.data.frame(SNP = Y$rsid, beta = Y$beta, se = Y$se, effect_allele = Y$A1, other_allele = Y$A2, chr=Y$chr, Phenotype=OUT, tstat=Y$tstat )

  mr_ind=unlist(prune_X(mr_dataX$tstat,pval))
  print(length(mr_ind))
  if(length(mr_ind)==0){
    mr_ind=unlist(prune_X(mr_dataX$tstat,pval1))
    print(length(mr_ind))
  }
  mr_dataX = mr_dataX[mr_ind,]
  mr_dataY = mr_dataY[mr_ind,]
  
  #filter( ( beta.exposure - beta.outcome ) / sqrt( se.exposure^2 + se.outcome^2 ) > reverse_t_threshold ) %>%
  ind_keep=which((abs(mr_dataX$beta)-abs(mr_dataY$beta))/sqrt(mr_dataX$se^2+mr_dataY$se^2) > reverse_t_threshold)
  print(length(ind_keep))
  mr_dataX = mr_dataX[ind_keep,]
  mr_dataY = mr_dataY[ind_keep,]
  
  exp_dat <- format_data(mr_dataX, type="exposure")  ##same rows as mr_dataX
  clump_bin = snp_bin(mr_dataX,50000)
  
  #exp_dat2 <- clump_data(exp_dat)
  exp_data = c()
  for (x in 1:length(clump_bin)) {
    temp = exp_dat[exp_dat$SNP %in% clump_bin[[x]]$SNP,]
    temp1 = clump_data(temp)
    exp_data=rbind(exp_data,temp1)
  }
  
  dups=which(duplicated(exp_data$SNP)==TRUE)
  if(length(dups)>0){
    exp_dat2 = exp_data[-dups,]
  }else{
    exp_dat2 = exp_data
  }
  
  out_dat <- format_data(mr_dataY, type="outcome")
  out_dat2=out_dat[out_dat$SNP %in% exp_dat2$SNP,]
  
  exp_dat2=exp_dat2[order(exp_dat2$SNP),] 
  out_dat2=out_dat2[order(out_dat2$SNP),] 
  
  if(all(exp_dat2$`effect_allele.exposure` == out_dat2$`effect_allele.outcome`)){
    print("action=1")
    action = 1
  } else {
    print("action=2/3")
    aligned = which(exp_dat2$`effect_allele.exposure` == out_dat2$`effect_allele.outcome` &
                      exp_dat2$`other_allele.exposure` == out_dat2$`other_allele.outcome`)
    swapped = which(exp_dat2$`effect_allele.exposure` == out_dat2$`other_allele.outcome` &
                      exp_dat2$`other_allele.exposure` == out_dat2$`effect_allele.outcome`)
    exp_dat2[swapped,'beta.exposure']=exp_dat2[swapped,'beta.exposure']*-1
    exp_dat2 = exp_dat2[c(aligned,swapped),]
    out_dat2 = out_dat2[c(aligned,swapped),]
    action = 1  ## made sure all strands are okay
  }
  # Harmonise the exposure and outcome data
  dat <- harmonise_data(
    exposure_dat = exp_dat2, 
    outcome_dat = out_dat2, action = action
  )
  
  #Sensitivity - Q-test
  het <- mr_heterogeneity(dat)
  het$I2 = ((het$Q-het$Q_df)/het$Q)*100
  plei <- mr_pleiotropy_test(dat)
  
  # Perform MR
  smaller=FALSE
  tryCatch( {res1 <- mr(dat); print("Bigger MR list") }, error = function(e) {smaller<<-TRUE})
  if(smaller){
    print("Smaller MR list")
    res1 <- mr(dat, method_list=c("mr_egger_regression","mr_weighted_median","mr_ivw")) 
  }else{
    res1 <- mr(dat)
  }
  
  alp_MR = res1[which(res1$method=="Inverse variance weighted"),'b']
  write.table(as.data.frame(res1), file=file_output, sep = ",", append = TRUE, row.names = FALSE)
  write.table(as.data.frame(het), file=file_output, sep = ",", append = TRUE, row.names = FALSE)
  write.table(as.data.frame(plei), file=file_output, sep = ",", append = TRUE, row.names = FALSE)
  write.table("*", file_output, sep = ",", quote = FALSE, col.names = FALSE, row.names = FALSE, append = TRUE)

  ## reverse MR
  rm(mr_dataX,mr_dataY,exp_dat,exp_data,exp_dat2,out_dat,out_dat2,dups,ind_keep,mr_ind,clump_bin, action, temp, temp1,dat,het,plei,smaller)
  
  #reverse the exposure and outcome to Y - X, nothing else besides this needs to change
  mr_dataX = cbind.data.frame(SNP = Y$rsid, beta = Y$beta, se = Y$se, effect_allele = Y$A1, other_allele = Y$A2, chr=Y$chr, Phenotype=OUT, tstat = Y$tstat )
  mr_dataY = cbind.data.frame(SNP = X$rsid, beta = X$beta, se = X$se, effect_allele = X$A1, other_allele = X$A2, chr=X$chr, Phenotype=EXP, tstat = X$tstat )
  
  mr_ind=unlist(prune_X(mr_dataX$tstat,pval))
  print(length(mr_ind))
  if(length(mr_ind)==0){
    mr_ind=unlist(prune_X(mr_dataY$tstat,pval1))
    print(length(mr_ind))
  }
  mr_dataX = mr_dataX[mr_ind,]
  mr_dataY = mr_dataY[mr_ind,]
  
  #filter( ( beta.exposure - beta.outcome ) / sqrt( se.exposure^2 + se.outcome^2 ) > reverse_t_threshold ) %>%  ind_keep=which((mr_dataX$beta-mr_dataY$beta)/sqrt(mr_dataX$se^2+mr_dataY$se^2) > reverse_t_threshold)
  #ind_keep=which((mr_dataX$beta-mr_dataY$beta)/sqrt(mr_dataX$se^2+mr_dataY$se^2) > reverse_t_threshold)
  ind_keep=which((abs(mr_dataX$beta)-abs(mr_dataY$beta))/sqrt(mr_dataX$se^2+mr_dataY$se^2) > reverse_t_threshold)
  print(length(ind_keep))
  mr_dataX = mr_dataX[ind_keep,]
  mr_dataY = mr_dataY[ind_keep,]
  
  exp_dat <- format_data(mr_dataX, type="exposure")  ##same rows as mr_dataX
  clump_bin = snp_bin(mr_dataX,50000)
  
  #exp_dat2 <- clump_data(exp_dat)
  exp_data = c()
  for (x in 1:length(clump_bin)) {
    temp = exp_dat[exp_dat$SNP %in% clump_bin[[x]]$SNP,]
    temp1 = clump_data(temp)
    exp_data=rbind(exp_data,temp1)
  }
  
  dups=which(duplicated(exp_data$SNP)==TRUE)
  if(length(dups)>0){
    exp_dat2 = exp_data[-dups,]
  }else{
    exp_dat2 = exp_data
  }
  
  out_dat <- format_data(mr_dataY, type="outcome")
  out_dat2=out_dat[out_dat$SNP %in% exp_dat2$SNP,]
  
  exp_dat2=exp_dat2[order(exp_dat2$SNP),] 
  out_dat2=out_dat2[order(out_dat2$SNP),] 
  
  if(all(exp_dat2$`effect_allele.exposure` == out_dat2$`effect_allele.outcome`)){
    print("action=1")
    action = 1
  } else {
    print("action=2/3")
    aligned = which(exp_dat2$`effect_allele.exposure` == out_dat2$`effect_allele.outcome` &
                      exp_dat2$`other_allele.exposure` == out_dat2$`other_allele.outcome`)
    swapped = which(exp_dat2$`effect_allele.exposure` == out_dat2$`other_allele.outcome` &
                      exp_dat2$`other_allele.exposure` == out_dat2$`effect_allele.outcome`)
    exp_dat2[swapped,'beta.exposure']=exp_dat2[swapped,'beta.exposure']*-1
    exp_dat2 = exp_dat2[c(aligned,swapped),]
    out_dat2 = out_dat2[c(aligned,swapped),]
    action = 1  ## made sure all strands are okay
  }
  
  # Harmonise the exposure and outcome data
  dat <- harmonise_data(
    exposure_dat = exp_dat2, 
    outcome_dat = out_dat2, action = action
  )
  
  #Sensitivity - Q-test
  het <- mr_heterogeneity(dat)
  het$I2 = ((het$Q-het$Q_df)/het$Q)*100
  plei <- mr_pleiotropy_test(dat)
  
  # Perform MR
  smaller=FALSE
  tryCatch( {res2 <- mr(dat); print("Bigger MR list") }, error = function(e) {smaller<<-TRUE})
  if(smaller){
    print("Smaller MR list")
    res2 <- mr(dat, method_list=c("mr_egger_regression","mr_weighted_median","mr_ivw")) 
  }else{
    res2 <- mr(dat)
  }
  
  bet_MR = res2[which(res2$method=="Inverse variance weighted"),'b']
  write.table(as.data.frame(res2), file=file_output, sep = ",", append = TRUE, row.names = FALSE)
  write.table(as.data.frame(het), file=file_output, sep = ",", append = TRUE, row.names = FALSE)
  write.table(as.data.frame(plei), file=file_output, sep = ",", append = TRUE, row.names = FALSE)
  
  return(list(i_XY, alp_MR, bet_MR))
  
}

  
