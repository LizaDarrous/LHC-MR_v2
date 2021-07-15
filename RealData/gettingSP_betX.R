## Per trait, getting the polygencity, direct heratibility and trait intercept (iX)

#library(psych)
library(data.table)
library(rslurm)
library(stringr)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
Gdir <- as.character(args[1])
EXPdir <- as.character(args[2])
trait <- as.character(args[3])
#partition <- as.character(args[4]) ### Boolean, change to false if complete model was already ran.

#Gdir = "/data/sgg3/liza/SEM_Realv2/data"
#EXPdir = "/data/sgg3/liza/SEM_Realv2/data/BMI_uniq.tsv"
#trait = "BMI"
#partition = "cluster2"

### updated function from rslurm that checks if slurm jobs are still running on cluster
get_job_status <- function (slr_job) 
{
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

grand_dir=Gdir ### Directory specific to data
#Should contain the summary stat files
setwd(grand_dir)
Xfile = fread(EXPdir) ### Exposure file if non-UKBB (overlapping SNPs between EXP-OUT)

LDfile = fread("/data/sgg3/liza/SEM_Realv2/data/LDscores.txt", sep="\t", header=TRUE)
## Extract SNPs high-quality SNPs, defined as being present in both UK10K and UK Biobank, having MAF>1 in both data sets,
#info>0.99 in the UK Biobank, non-significant (P_{diff}>0.05) allele frequency difference between UK Biobank and UK10K and
#residing outside the HLA region (chr6:28.5-33.5Mb)
attach(LDfile)
mafT = .005
selF = which(info>.99 & abs(mafUK10K-mafUKBB)<(2*sqrt(1/4e3+1/360e3)) & mafUKBB>mafT & mafUK10K>mafT &
               !(chr==6 & pos>=28.5e6 & pos<=33.5e6))
LDfile = LDfile[selF,]

### get pi1 and sigma1 and SNPs in common
rho_file = fread("/data/sgg2/zoltan/project/Project_SEM/results/LD_GM2_2prm.txt", header=FALSE)
colnames(rho_file) = c("rsid", "chr", "pos", "piK", "sigK")

Xfile %>% arrange(chr, pos) -> X_data

Data = inner_join(X_data, rho_file) # based on rsid / chr / pos
nrow(Data)
# 4,689,924

colnames(LDfile)[3] = "rsid"
LDfile$info = NULL # lots of SNPs have different info (UKBB vs UK10K?), needed otherwise only 207,914 SNPs left
Data = inner_join(Data, LDfile) # based on rsid / chr / pos
nrow(Data)
# 4,689,924

# slice, every 10th
Data %>%
  slice(seq(1, nrow(Data), by=10)) -> Data_filtered
nrow(Data_filtered)
# 468,993

nX = mean(Data_filtered$n_complete_samples)  #Get sample size for trait X

bX = Data_filtered$tstat/sqrt(nX)   #Get standardised beta for trait X

ld = Data_filtered$LDSC
weights = Data_filtered$weight
pi1 = Data_filtered$piK
sig1 = Data_filtered$sigK

betX = bX #cbind(bX, bY)  # used directly in LHC_SEM as a global parameter
# BE CAREFUL, then, if you modify it in LHC_SEM, it will be definetively changed???
m0 = (2500*2)+1
M = 10106833  #Number of SNPs

sp_piX = runif(50,0,0.01)
sp_h2X = runif(50,0,0.5)
sp_iX = runif(50,0.5,1.5)

para=cbind(sp_piX,sp_h2X,sp_iX)
sp_mat=matrix(unlist(para), ncol=3, byrow = FALSE)
colnames(sp_mat)=colnames(para)

### likelihood function
run_optim_noHess = function(par){
  source("/data/sgg3/liza/SEM_Realv2/scripts/LHC_MR_betX_v9_LD.R")
  
  #theta= sp_mat[1,]
  #LHC_MR_betX_v9_LD(theta,betX,pi1,sig1,weights,m0,nX,bn=2^7,bins=10)
  
  theta=unlist(par)
  test = optim(theta, LHC_MR_betX_v9_LD,
               betX=betX, pi1=pi1, sig1=sig1, weights=weights,
               m0=m0, nX=nX, bn=2^7, bins=10,
               method = "Nelder-Mead",
               #lower = c(0,0,0,0,0,-1,-1,-1,0,0,-1),
               #upper = c(1,1,1,1,1,1,1,1,2,2,1),
               control = list(maxit = 5e3))
  
  test.res=c(test$value,test$par,test$convergence)
  names(test.res)=c("mLL","piX", "h2X","iX","conv")
  return(test.res)
}

### optimisation
par.df = data.frame(par=I(apply(sp_mat,1,as.list))) #generate a dataframe of lists for each row of parameters - input for rslurm

start.time <- Sys.time()
sjob = slurm_apply(f = run_optim_noHess, params = par.df, jobname = paste0("SP_",trait), nodes = 100, cpus_per_node = 1,
                   #libPaths = c(.libPaths(), "/data/sgg2/ninon/bin/R-3.4.3_Packages/"),
                   add_objects = c("betX","pi1","sig1","weights","m0","nX"),
                   slurm_options = list(partition = "sgg"), # , `cpus-per-task`=4 #note: partition X is representative of a single partition in cluster.
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
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)

value.mat = get_slurm_out(sjob, outtype = 'table')
#res = cbind(sp_mat,value.mat)
#write.table(res, file = paste0(trait,"_SP_detailed.csv"), sep = ",", quote = FALSE, col.names = TRUE, row.names = FALSE)

res_min = value.mat[which(value.mat$mLL == min(value.mat$mLL)), ]
res = cbind(trait,abs(res_min[2:4]))
write.table(res, file = "SingleTrait_SP.csv", sep = ",", quote = FALSE, col.names = FALSE, row.names = FALSE, append = TRUE)
write.table("*", file = "SingleTrait_SP.csv", sep = ",", quote = FALSE, col.names = FALSE, row.names = FALSE, append = TRUE)

cleanup_files(sjob)
