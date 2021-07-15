#!/bin/bash

#SBATCH --chdir=/project/scripts # The Working Directory of the job

Rscript DataOptimisation_Automated_7params.R /project/results /project/data EXPOSURE_name OUTCOME_name PARTITION
