#!/bin/bash

#SBATCH --chdir=/data/sgg3/liza/SEM_Realv2/scripts # The Working Directory of the job

Rscript DataOptimisation_Automated_7params.R /data/sgg3/liza/SEM_Realv2/resultsv2 /data/sgg3/liza/SEM_Realv2/data HDL Asthma sgg
Rscript DataOptimisation_Automated_7params.R /data/sgg3/liza/SEM_Realv2/resultsv5 /data/sgg3/liza/SEM_Realv2/data HDL DM sgg
Rscript DataOptimisation_Automated_7params.R /data/sgg3/liza/SEM_Realv2/resultsv8 /data/sgg3/liza/SEM_Realv2/data HDL SBP sgg
Rscript DataOptimisation_Automated_7params.R /data/sgg3/liza/SEM_Realv2/resultsv11 /data/sgg3/liza/SEM_Realv2/data HDL MI sgg
Rscript DataOptimisation_Automated_7params.R /data/sgg3/liza/SEM_Realv2/resultsv14 /data/sgg3/liza/SEM_Realv2/data LDL BMI sgg
Rscript DataOptimisation_Automated_7params.R /data/sgg3/liza/SEM_Realv2/resultsv17 /data/sgg3/liza/SEM_Realv2/data LDL Edu sgg
Rscript DataOptimisation_Automated_7params.R /data/sgg3/liza/SEM_Realv2/resultsv20 /data/sgg3/liza/SEM_Realv2/data LDL SHeight sgg
Rscript DataOptimisation_Automated_7params.R /data/sgg3/liza/SEM_Realv2/resultsv23 /data/sgg3/liza/SEM_Realv2/data CAD Asthma sgg
Rscript DataOptimisation_Automated_7params.R /data/sgg3/liza/SEM_Realv2/resultsv26 /data/sgg3/liza/SEM_Realv2/data CAD DM sgg
Rscript DataOptimisation_Automated_7params.R /data/sgg3/liza/SEM_Realv2/resultsv29 /data/sgg3/liza/SEM_Realv2/data CAD SBP sgg
Rscript DataOptimisation_Automated_7params.R /data/sgg3/liza/SEM_Realv2/resultsv32 /data/sgg3/liza/SEM_Realv2/data CAD MI sgg
Rscript DataOptimisation_Automated_7params.R /data/sgg3/liza/SEM_Realv2/resultsv35 /data/sgg3/liza/SEM_Realv2/data DM Asthma sgg
Rscript DataOptimisation_Automated_7params.R /data/sgg3/liza/SEM_Realv2/resultsv38 /data/sgg3/liza/SEM_Realv2/data DM SBP sgg
Rscript DataOptimisation_Automated_7params.R /data/sgg3/liza/SEM_Realv2/resultsv41 /data/sgg3/liza/SEM_Realv2/data DM MI sgg
Rscript DataOptimisation_Automated_7params.R /data/sgg3/liza/SEM_Realv2/resultsv44 /data/sgg3/liza/SEM_Realv2/data SHeight MI sgg
Rscript DataOptimisation_Automated_7params.R /data/sgg3/liza/SEM_Realv2/resultsv47 /data/sgg3/liza/SEM_Realv2/data BMI DM sgg
Rscript DataOptimisation_Automated_7params.R /data/sgg3/liza/SEM_Realv2/resultsv50 /data/sgg3/liza/SEM_Realv2/data BMI SBP sgg
Rscript DataOptimisation_Automated_7params.R /data/sgg3/liza/SEM_Realv2/resultsv53 /data/sgg3/liza/SEM_Realv2/data BMI MI sgg
Rscript DataOptimisation_Automated_7params.R /data/sgg3/liza/SEM_Realv2/resultsv56 /data/sgg3/liza/SEM_Realv2/data BWeight Edu sgg
Rscript DataOptimisation_Automated_7params.R /data/sgg3/liza/SEM_Realv2/resultsv59 /data/sgg3/liza/SEM_Realv2/data BWeight SHeight sgg
Rscript DataOptimisation_Automated_7params.R /data/sgg3/liza/SEM_Realv2/resultsv62 /data/sgg3/liza/SEM_Realv2/data Edu Asthma sgg
Rscript DataOptimisation_Automated_7params.R /data/sgg3/liza/SEM_Realv2/resultsv65 /data/sgg3/liza/SEM_Realv2/data Edu SHeight sgg
Rscript DataOptimisation_Automated_7params.R /data/sgg3/liza/SEM_Realv2/resultsv68 /data/sgg3/liza/SEM_Realv2/data PSmoke Asthma sgg
Rscript DataOptimisation_Automated_7params.R /data/sgg3/liza/SEM_Realv2/resultsv71 /data/sgg3/liza/SEM_Realv2/data PSmoke SVstat sgg
Rscript DataOptimisation_Automated_7params.R /data/sgg3/liza/SEM_Realv2/resultsv77 /data/sgg3/liza/SEM_Realv2/data SVstat Asthma sgg