#!/bin/bash

#SBATCH --chdir=/project/scripts # The Working Directory of the job


Rscript gettingSP_betX.R /project/data/ /project/data/TRAIT_uniq.tsv TRAIT_name

