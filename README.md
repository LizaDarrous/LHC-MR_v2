# LHC-MR_v2
(Updated method)

:warning: **These scripts are in the process of being made into an R-package. Thank you for your patience!**

Latent Heritable Confounder MR (LHC-MR) is a method applicable to association summary statistics, which estimates bi-directional causal effects, direct heritabilities, and confounder effects while accounting for sample overlap.
LHC-MR extends the standard Mendelian Randomisation model to incorporate the presence of a latent (unmeasured) heritable confounder and estimates its contribution to the exposure and outcome traits, while simultaneously estimating the bi-directional causal effect between the two traits.


## Usage

The source code for LHC-MR was written in R version 4.0.2. No other language is used in the computation, and thus only bash and R are needed. 
Several R-packages are needed to run the analysis, which will be detailed below.


### Association Summary Statistics
#### Prerequisites
R Packages needed to run the summary stat analysis include:
``` 
install.packages(data.table);  library(data.table)
install.packages(DescTools); library(DescTools)
install.packages(rslurm);  library(rslurm)
install.packages(tidyverse); library(tidyverse)
install.packages(stringr); library(stringr)
install.packages(tictoc);  library(tictoc)
install.packages(psych);  library(psych)
install.packages(extRemes);  library(extRemes)
install.packages("remotes");remotes::install_github("MRCIEU/TwoSampleMR");library(TwoSampleMR)
remotes::install_github("GenomicSEM/GenomicSEM");library(GenomicSEM)

install.packages(GGally);  library(GGally)
install.packages(mixtools);  library(mixtools)
```
A very important package used in the analysis is [rslurm](https://cran.r-project.org/web/packages/rslurm/rslurm.pdf), which allows us to submit array jobs from within R without having to create a bash script to do so. Moreover, this parallelisation step takes advantage of the presence of a cluster with several partitions, onto which it can simultaneously distribute and run array jobs. Once rslurm is installed, it's important to edit the template files that come with the package to reflect the info of the cluster being used.

The LHC-MR scripts was written with a basic folder hierarchy assummed, such as the following (note that this is just a **suggestion**), but the paths can be changed in the scripts currently to better reflect the users folder hierarchy:
```
project --|---- data
          |       |--- Raw/Unprocessed GWAS summary stats files
          |       |--- Processed GWAS summary stats files
          |       |--- variants.tsv.bgz
          |       |--- LDscores.txt
          |       |--- LD_GM_2prm.txt                   
          |       |---- genomicSEM (folder)
          |
          |---- scripts
          |
          |---- results

```
Raw summary stat files downloaded in `data` require certain columns to continue with analysis: rsid, EffectAllele/A1/alt, OtherAllele/A2/ref, number of samples (effective sample size), beta/effect size of effective allele,
se/standard error of effect size.

These files ought to be processed using the `get_tX_UKBB_Automated.R` script run through bash for bulk traits, in order to calculate the t-statistics (beta/se) of the SNPs, which will be used in the main LHC-MR analysis. This code can be generalised for traits coming from different cohorts.

This script first selects the set of SNPs for which imputation quality info>0.99 and MAF>0.005 using the `LDscore.txt` file. 
For LD scores you can either use the classical ones (https://data.broadinstitute.org/alkesgroup/LDSCORE/) or we have also created LD scores based on UK10K sequencing data (https://drive.google.com/file/d/1uua6zIournPvcJAs-QVWs8pTUVLtZ5Ym/view?usp=sharing), where the list of SNPs are optimised for UKB-based summary stats (e.g. by Neale).
Then for these set of SNPs, it uses the [`variants.tsv.bgz`](https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/variants.tsv.bgz) obtained from Neale's [UKBB GWAS Imputed v3](https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit?ts=5b5f17db#gid=178908679) to merge variant information to variant details of the summary statistics files. This leads to the addition of the A1 and A2 columns as well as the rsid, all of which are needed for further analysis. Some duplication occurs when `variants.tsv.bgz` is merged with UKBB summary stat files, which is handled in the script.
This script is used once per UKBB trait, allowing all the traits to then have the same set of overlapping SNPs when running trait pair analysis between UKBB traits. 

Following this step, a single trait LHC-MR analysis is ran using the script `gettingSP_betX.R` to obtain values for the poligencity and single trait intercept, from which starting points will be derived in the pair trait analysis. This can be run for all traits you wish to analyse and will store these values in a file written in the `data` folder which will be read in `DataOptimisation_Automated_7params.R` when needed.
The parameters estimated here are:-
- iX: single trait intercept, representative of the population structure.
- pX: poligencity of X (proportion of SNPs with an effect on trait X).
- h2X: total heritabilities of X.

For the LHC-MR analysis, the script `DataOptimisation_Automated_7params.R` first reads in the summary stat files for the exposure (X) and outcome (Y). It then has the option to run LDSC using the `genomicSEM` R package and standardMR methods using the `TwoSampleMR` R package in order to obtain a cross trait intercept and causal effects that can be used as starting points for the likleihood optimisation. Random number sare generated if this step is skipped.
Specific files are needed to run `genomicSEM`, and more details on abaining them can be found in the `data` folder above.
The likelihood function is estimated with inputs including the sample sizes, the standardised effect sizes, the SNP-specific LD information (obtained from `LD_GM2_2prm.txt`), and the starting points.
7 unknown parameters are then estimated:-
- tX, tY: confounder effect of X and Y.
- alp, bet: causal effect of X->Y and Y->X.
- rho: representative of the phenotypic correlation due to sample overlap.
The results are reported in the trait-subfolder that's created called `TwoStep`

Once these 7 parameters are estimated, a 200 block jacknife step follows to estimate the SE of these parameters. The detailed results will be created in the trait-subfolder `block200JK1sp`.

## Citation

Preprint: https://www.medrxiv.org/content/10.1101/2020.01.27.20018929v2.article-metrics 

## Authors

Liza Darrous - liza.darrous@unil.ch

Ninon Mounier

Zoltán Kutalik

## Acknowledgments

SGG Lab members

Package creators for those used
