#!/bin/bash

#SBATCH --chdir=/project/scripts # The Working Directory of the job

Rscript get_tX_UKBB_Automated.R /project/data/20002_1111.gwas.imputed_v3.both_sexes.tsv.bgz TRAIT_name
