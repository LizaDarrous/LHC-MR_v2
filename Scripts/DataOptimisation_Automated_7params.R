### Two step optimisation - smaller bins first then using the top few results as SP for bigger bins
library(tictoc)
library(psych)
library(data.table)
library(rslurm)
library(stringr)
library(tidyverse)
library(extRemes)

args <- commandArgs(trailingOnly = TRUE)
Rdir <- as.character(args[1]) #result folder "/project/results"
Ddir <- as.character(args[2]) #data folder "/project/data"
EXP <- as.character(args[3])
OUT <- as.character(args[4])
partition <- as.character(args[5]) #if the cluster has several partitions

### updated function from rslurm that checks if slurm jobs are still running on cluster
get_job_status <- function (slr_job){
  if (!(class(slr_job) == "slurm_job")) 
    stop("input must be a slurm_job")
  stat <- suppressWarnings(system(paste("squeue -n", slr_job$jobname), 
                                  intern = TRUE))
  if (length(stat) > 1) {
    res = "Job running or in queue."
  }
  else {
    res = "Job completed or stopped."
    # tmpdir <- paste0("rslurm", slr_job$jobname)
    #  out_files <- file.path(tmpdir, paste0("slurm_", 0:(slr_job$nodes - 
    #                                                      1), ".out"))
    #for (outf in out_files) {
    #  cat(paste("\n----", outf, "----\n\n"))
    #  cat(paste(readLines(outf), collapse = "\n"))
  }
  return(res)
}

### likelihood function
run_optim_Hess = function(par){
  source("/project/scripts/LHC_MR_bXY_v91_LD.R") ## likelihood function
  theta=unlist(par)
  
  test = optim(theta, LHC_MR_bXY_v91_LD,
               betXY=betXY, pi1=pi1, sig1=sig1, weights=weights, pi_U=pi_U,
               pi_X=pi_X, pi_Y=pi_Y, i_X=i_X, i_Y=i_Y,
               m0=m0, nX=nX, nY=nY, bn=bn, bins=bins, model=param,
               method = "Nelder-Mead",
               hessian = TRUE,
               #lower = c(0,0,0,0,0,-1,-1,-1,0,0,-1),
               #upper = c(1,1,1,1,1,1,1,1,2,2,1),
               control = list(maxit = 5e3,
                              parscale = parscale))
  
  test.res=c(test$value,test$par,test$convergence)
  test.se = c(test$value,sqrt(diag(solve(test$hessian))))
  cnames_res = c("mLL","h2X","h2Y","tX","tY","alp","bet","iXY","conv")
  cnames_se = c("mLL","h2X","h2Y","tX","tY","alp","bet","iXY")
  
  if(param=="comp"){
    names(test.res)=cnames_res
    names(test.se)=cnames_se
  }else if(param=="U"){
    param_temp = c("tX","tY")
    names(test.res)=cnames_res[!(cnames_res %in% param_temp)]
    names(test.se)=cnames_se[!(cnames_se %in% param_temp)]
  }else{
    names(test.res)=cnames_res[!(cnames_res %in% param)]
    names(test.se)=cnames_se[!(cnames_se %in% param)]
  }
  
  return(list(res = test.res, se = test.se))
}

setwd(Rdir)
pair_dir = paste0(EXP,"-",OUT)
system(paste0("mkdir ./",pair_dir))
setwd(pair_dir)
system("mkdir ./TwoStep")
grand_dir=paste0(Rdir,"/",pair_dir,"/TwoStep") ### Directory specific to the trait pair to be analysed.
#Should contain the summary stat files if non-UKBB SNPs were made to overlap in this directory.
setwd(grand_dir)

Xfile = fread(paste0(Ddir,"/",EXP,"_uniq.tsv")) #fread(EXPdir) ### Exposure file if non-UKBB (overlapping SNPs between EXP-OUT)
Yfile = fread(paste0(Ddir,"/",OUT,"_uniq.tsv")) #fread(OUTdir) ### Outcome file if non-UKBB (overlapping SNPs between EXP-OUT)

LDfile = fread("/project/data/LDscores.txt", sep="\t", header=TRUE)
## Extract SNPs high-quality SNPs, defined as being present in both UK10K and UK Biobank, having MAF>1 in both data sets, 
#info>0.99 in the UK Biobank, non-significant (P_{diff}>0.05) allele frequency difference between UK Biobank and UK10K and 
#residing outside the HLA region (chr6:28.5-33.5Mb)
attach(LDfile)
mafT = .005
selF = which(info>.99 & abs(mafUK10K-mafUKBB)<(2*sqrt(1/4e3+1/360e3)) & mafUKBB>mafT & mafUK10K>mafT & 
               !(chr==6 & pos>=28.5e6 & pos<=33.5e6))
LDfile = LDfile[selF,]

### get pi1 and sigma1 and SNPs in common
rho_file = fread("/project/data/LD_GM2_2prm.txt", header=FALSE)
colnames(rho_file) = c("rsid", "chr", "pos", "piK", "sigK")

# slightly change pre-processing code + order by chr/pos before slicing!
Xfile %>% arrange(chr, pos) -> X_data
Yfile %>% arrange(chr, pos) -> Y_data
Data = inner_join(X_data, Y_data,
                  by = c("chr", "pos", "rsid"))

aligned = which(Data$alt.x==Data$alt.y &
                  Data$ref.x==Data$ref.y)
swapped = which(Data$alt.x==Data$ref.y &
                  Data$ref.x==Data$alt.y)
#Correct the effect of swapped alleles as well as the t-stat for one of the two traits
Data[swapped,'tstat.x']=Data[swapped,'tstat.x']*-1
Data[swapped,'beta.x']=Data[swapped,'beta.x']*-1
temp_alt=Data$alt.x
Data[swapped,'alt.x']=Data[swapped,'ref.x']
Data[swapped,'ref.x']=temp_alt[swapped]

Data1=Data[c(aligned,swapped),]
## test swapping
all(Data1$alt.x==Data1$alt.y)
all(Data1$ref.x==Data1$ref.y)
Data = Data1
Data$alt = Data$alt.x
Data$ref = Data$ref.x
nrow(Data) # 4,732,967

Data = inner_join(Data, rho_file) # based on rsid / chr / pos
nrow(Data) # 4,689,922

colnames(LDfile)[3] = "rsid"
LDfile$info = NULL # lots of SNPs have different info (UKBB vs UK10K?), needed otherwise only 207,914 SNPs left
Data = inner_join(Data, LDfile) # based on rsid / chr / pos
nrow(Data) # 4,689,922

# slice, every 10th
Data %>%
  slice(seq(1, nrow(Data), by=10)) -> Data_filtered
nrow(Data_filtered) # 468,993

nX = mean(Data_filtered$n_complete_samples.x)  #Get sample size for trait X
nY = mean(Data_filtered$n_complete_samples.y)  #Get sample size for trait Y

bX = Data_filtered$tstat.x/sqrt(nX)   #Get standardised beta for trait X
bY = Data_filtered$tstat.y/sqrt(nY)   #Get standardised beta for trait Y

ld = Data_filtered$LDSC
weights = Data_filtered$weight
pi1 = Data_filtered$piK
sig1 = Data_filtered$sigK

betXY = cbind(bX, bY)  # used directly in LHC_SEM as a global parameter
m0 = (2500*2)+1
#M = 10106833  #Number of SNPs

### get starting points + generate the rest
betX_df=fread("/project/data/SingleTrait_SP.csv") ## Estimates of pi, h2, and iX from SingleTrait analysis (gettingSP_betX.R)
colnames(betX_df)=c("Trait","piX","h2X","iX")
i_X = unlist(betX_df[betX_df$Trait==EXP,"iX"])
i_Y = unlist(betX_df[betX_df$Trait==OUT,"iX"])
h2_x = unlist(betX_df[betX_df$Trait==EXP,"h2X"])
h2_y = unlist(betX_df[betX_df$Trait==OUT,"h2X"])
pi_X = unlist(betX_df[betX_df$Trait==EXP,"piX"])
pi_Y = unlist(betX_df[betX_df$Trait==OUT,"piX"])

source("/project/scripts/gettingSP_ldscMR.R")

i_XY = SP[[1]]
alp_MR = SP[[2]]
bet_MR = SP[[3]]

sp_tX = runif(100,0,0.5)
sp_tY = runif(100,-0.5,0.5)
sp_h2X = h2_x-(sp_tX^2)
sp_h2Y = h2_y-(sp_tY^2)
sp_alp = replicate(100, (alp_MR+runif(1,-0.1,0.1)))
sp_bet = replicate(100, (bet_MR+runif(1,-0.1,0.1)))
sp_iXY = rep(i_XY,100)

para=cbind(sp_h2X,sp_h2Y,sp_tX,sp_tY,sp_alp,sp_bet,sp_iXY)
sp_mat=matrix(unlist(para), ncol=7, byrow = FALSE)
colnames(sp_mat)=colnames(para)
sp_mat1 = cbind("SP"=c(1:nrow(sp_mat)),sp_mat)
write.csv(sp_mat1,"100_startingpoints_Small.csv", row.names=F)

### Small grid optimisation
par.df = data.frame(par=I(apply(sp_mat,1,as.list))) #generate a dataframe of lists for each row of parameters - input for rslurm
pi_U=0.1
bn = 2^7
bins = 10
param="comp"
parscale = c(1e-1,1e-1,1e-1,1e-1,1e-1,1e-1,1e-1)
tic()
sjob = slurm_apply(f = run_optim_Hess, params = par.df, jobname = paste0(EXP,OUT,"-Small"), nodes = 100, cpus_per_node = 1,
                   #libPaths = c(.libPaths(), "/data/sgg2/ninon/bin/R-3.4.3_Packages/"),
                   add_objects = c("betXY","pi1","sig1","weights","m0","nX","nY","pi_U","pi_X","pi_Y",
                                   "i_X","i_Y","param","bn","bins","parscale"),
                   slurm_options = list(partition = partition, time="1-00:00:00"),
                   submit = TRUE)
#Keep a loop open till the job is completely done.
wait_counter = 0
while (wait_counter < 1) {
  wait_counter = 0
  if (tryCatch(str_detect(get_job_status(sjob),"completed"), warning = function(w){FALSE},
               error = function(e){FALSE})) {
    wait_counter = wait_counter + 1
  } else{
    wait_counter = wait_counter
  }
}
time.taken = toc() #4253.686 sec elapsed
print(time.taken)

res_temp = get_slurm_out(sjob, outtype = 'table')
res_values = as.data.frame(do.call(rbind, res_temp[[1]]))
res_values %>%
  mutate(h2X = abs(h2X),
         h2Y = abs(h2Y),
         tX = abs(tX)) -> res_values
res_values = cbind("SP"=c(1:nrow(res_values)),"mLL"=res_values[,1],"piX"=pi_X,"piY"=pi_Y,res_values[,-1],"iX"=i_X,"iY"=i_Y)

res_se = as.data.frame(do.call(rbind, res_temp[[2]]))
res_se = cbind("SP"=c(1:nrow(res_se)),res_se)

write.csv(res_values, "AllRes_Hess_Small.csv", row.names = FALSE) ## Will be used to read new sp_mat
write.csv(res_se, "AllRes_Hess_SE_Small.csv", row.names = FALSE) ## Will be used to read new sp_mat

cleanup_files(sjob)

#### block JK
### likelihood function
run_optim_Hess_JK = function(par,start_ind, end_ind){
  source("/project/scripts/LHC_MR_bXY_v91_LD.R")
   
  theta=unlist(par)
  
  test = optim(theta, LHC_MR_bXY_v91_LD,
               betXY=betXY[-(start_ind:end_ind),], pi1=pi1[-(start_ind:end_ind)], sig1=sig1[-(start_ind:end_ind)],
               weights=weights[-(start_ind:end_ind)], pi_U=pi_U,
               pi_X=pi_X, pi_Y=pi_Y, i_X=i_X, i_Y=i_Y,
               m0=m0, nX=nX, nY=nY, bn=bn, bins=bins, model=param,
               method = "Nelder-Mead",
               hessian = TRUE,
               #lower = c(0,0,0,0,0,-1,-1,-1,0,0,-1),
               #upper = c(1,1,1,1,1,1,1,1,2,2,1),
               control = list(maxit = 5e3,
                              parscale = parscale))
  
  test.res=c(test$value,test$par,test$convergence,start_ind, end_ind)
  test.se = c(test$value,sqrt(diag(solve(test$hessian))),start_ind, end_ind)
  cnames_res = c("mLL","h2X","h2Y","tX","tY","alp","bet","iXY","conv","start_ind", "end_ind")
  cnames_se = c("mLL","h2X","h2Y","tX","tY","alp","bet","iXY","start_ind", "end_ind")
  
  if(param=="comp"){
    names(test.res)=cnames_res
    names(test.se)=cnames_se
  }else if(param=="U"){
    param_temp = c("tX","tY")
    names(test.res)=cnames_res[!(cnames_res %in% param_temp)]
    names(test.se)=cnames_se[!(cnames_se %in% param_temp)]
  }else{
    names(test.res)=cnames_res[!(cnames_res %in% param)]
    names(test.se)=cnames_se[!(cnames_se %in% param)]
  }
  
  return(list(res = test.res, se = test.se))
}

setwd(Rdir)
pair_dir = paste0(EXP,"-",OUT)
setwd(pair_dir)
system("mkdir ./block200JK1sp")
grand_dir=paste0(Rdir,"/",pair_dir,"/block200JK1sp") ### Directory specific to the trait pair to be analysed.
setwd(grand_dir)

sink(paste0(EXP,"-",OUT,"_log.txt"), append=FALSE, split=TRUE)

### getting blocks
uniqPos = seq(1:nrow(Data_filtered))
### Divide the snps into 202 blocks with 11 snps leftover added to the last chunk
nBlock = 200
nSNP = nrow(betXY) %/% nBlock  #2344
limit = nBlock * nSNP #468800
leftover = as.numeric(nrow(betXY) %% nSNP ) #193
start_ind = uniqPos[seq(1, limit, nSNP)]
end_ind = start_ind+(nSNP-1)
end_ind[length(end_ind)]=end_ind[length(end_ind)]+leftover

JK_index=cbind(start_ind,end_ind) # 200 rows/bins

res_ordered = res_values[order(res_values$mLL, decreasing = F),]
sp_mat = select(res_ordered,h2X,h2Y,tX,tY,alp,bet,iXY)
sp_mat = sp_mat[1,]

par.df = data.frame(par=I(apply(sp_mat,1,as.list))) #generate a dataframe of lists for each row of parameters - input for rslurm
par.df2 = merge(par.df,JK_index)
pi_U=0.1
bn = 2^7
bins = 10
param="comp"
parscale = c(1e-1,1e-1,1e-1,1e-1,1e-1,1e-1,1e-1)
tic()
sjob = slurm_apply(f = run_optim_Hess_JK, params = par.df2, jobname = paste0(EXP,OUT,"-blockJK"), nodes = 4000, cpus_per_node = 1,
                   add_objects = c("betXY","pi1","sig1","weights","m0","nX","nY","pi_U","pi_X","pi_Y",
                                   "i_X","i_Y","param","bn","bins","parscale"),
                   slurm_options = list(partition = partition,time="1-00:00:00", mem="3G"), # , `cpus-per-task`=4 #note: partition X is representative of a single partition in cluster.
                   submit = TRUE)
#Keep a loop open till the job is completely done.
wait_counter = 0
while (wait_counter < 1) {
  wait_counter = 0
  if (tryCatch(str_detect(get_job_status(sjob),"completed"), warning = function(w){FALSE},
               error = function(e){FALSE})) {
    wait_counter = wait_counter + 1
  } else{
    wait_counter = wait_counter
  }
}
time.taken = toc() #4253.686 sec elapsed
print(time.taken)

#Get output of minus log likelihood (mLL) and estimated parameters from rslurm in the form of a table with nrows equal to nrows(par.df)
res_temp = get_slurm_out(sjob, outtype = 'table')
res_values = as.data.frame(do.call(rbind, res_temp[[1]]))
#names(res_values) = c("mLL","piX","piU", "piY", "h2X","h2Y","tX","tY","alp","bet","iX","iY","iXY","conv")
res_values %>%
  mutate(h2X = abs(h2X),
         h2Y = abs(h2Y),
         tX = abs(tX)) -> res_values
#res_values = cbind("SP"=c(1:nrow(res_values)),res_values)
res_se = as.data.frame(do.call(rbind, res_temp[[2]]))
#names(res_se) = c("mLL","piX","piU", "piY", "h2X","h2Y","tX","tY","alp","bet","iX","iY","iXY")
#res_se = cbind("SP"=c(1:nrow(res_se)),res_se)

write.csv(res_values, "AllRes_Hess_JK.csv", row.names = FALSE) ## Will be used to read new sp_mat
write.csv(res_se, "AllRes_Hess_SE_JK.csv", row.names = FALSE) ## Will be used to read new sp_mat

res_min = res_values %>%
  group_by(start_ind) %>%
  slice(which.min(mLL)) %>%
  ungroup()

write.csv(res_min, "MinRes_JK.csv", row.names = FALSE) ## same as above with 1SP

res_minFil = dplyr::select(as.data.frame(res_min), -c(mLL,conv,start_ind,end_ind))
JK_res = as.data.frame(matrix(data=NA, nrow=ncol(res_minFil),ncol=18))
colnames(JK_res)=c("Parameter","mean","median","se","se_JK","se_JK_10","var","loglik_1comp","AIC_1comp","convergance_2comp","mu1","mu2","sigma1","sigma2","lambda1","lambda2","loglik_2comp","AIC_2comp")
JK_res$Parameter=colnames(res_minFil)
JK_res$mean = colMeans(res_minFil)
JK_res$median = apply(res_minFil, 2, median)
JK_res$se = apply(res_minFil, 2, sd)
JK_res$se_JK = (apply(res_minFil, 2, sd))*sqrt(nBlock-1)#*sqrt(10)
JK_res$se_JK_10 = (apply(res_minFil, 2, sd))*sqrt(nBlock-1)*sqrt(10)
JK_res$var = apply(res_minFil, 2, var)

for (x in c(1:ncol(res_minFil))){
  param = res_minFil[,x]
  print(colnames(res_minFil)[x])
  Xf1 = MASS::fitdistr(param, "normal")
  Xf2 = tryCatch({capture.output(mixtools::normalmixEM(param, k=2, maxit=1e8))}, warning = function(warning_condition) {
    return(NA)
  }, error = function(error_condition) {
    print("Error")
    return(NA)
  })
  AIC1 = 2*2 - 2*(Xf1$loglik)
  AIC2 = 2*4 - 2*(as.numeric(str_split(str_squish(Xf2[str_which(Xf2, "loglik")[1]+1]), " " )[[1]][2]))
  JK_res$loglik_1comp[x] = Xf1$loglik
  JK_res$loglik_2comp[x] = as.numeric(str_split(str_squish(Xf2[str_which(Xf2, "loglik")[1]+1]), " " )[[1]][2])
  JK_res$AIC_1comp[x] = AIC1
  JK_res$AIC_2comp[x] = AIC2
  
  
  if(!is.na(Xf2)){
    print(colnames(res_minFil)[x])
    JK_res$convergance_2comp[x] = !any(str_detect(Xf2, "WARNING! NOT CONVERGENT!"))
    JK_res[x,11:12] = as.numeric(str_split(str_squish(Xf2[str_which(Xf2, "mu")+1]), " " )[[1]][2:3])
    JK_res[x,13:14] = (as.numeric(str_split(str_squish(Xf2[str_which(Xf2, "sigma")+1]), " " )[[1]][2:3]))
    JK_res[x,15:16] = as.numeric(str_split(str_squish(Xf2[str_which(Xf2, "lambda")+1]), " " )[[1]][2:3])
  }
}

tstat = (JK_res$mu1 - JK_res$mu2) / sqrt(JK_res$sigma1^2 + JK_res$sigma2^2)
t.pval = 2*pnorm(-abs(tstat))

JK_res$tstat = tstat
JK_res$tstat_pval = t.pval
JK_res$sigma1_JK = JK_res$sigma1*sqrt(nBlock-1)
JK_res$sigma2_JK = JK_res$sigma2*sqrt(nBlock-1)
JK_res$sigma1_JK_10 = JK_res$sigma1*sqrt(nBlock-1)*sqrt(10)
JK_res$sigma2_JK_10 = JK_res$sigma2*sqrt(nBlock-1)*sqrt(10)

JK_res$ci_lower_JK = JK_res$mean - (1.96*JK_res$se_JK)
JK_res$ci_upper_JK = JK_res$mean + (1.96*JK_res$se_JK)
JK_res$ci_lower_JK_10 = JK_res$mean - (1.96*JK_res$se_JK_10)
JK_res$ci_upper_JK_10 = JK_res$mean + (1.96*JK_res$se_JK_10)

bimo = which(JK_res$tstat_pval == 0)

JK_res$bimod = "FALSE"
JK_res$bimod[bimo] = "TRUE"

write.csv(JK_res,paste0("JKres_200_",EXP,"-",OUT,".csv"), row.names = F)

cov_matrix = cov(res_minFil)
write.csv(cov_matrix,paste0("VarCovMatrix_",EXP,"-",OUT,".csv"))
print("Done!")
sink()
