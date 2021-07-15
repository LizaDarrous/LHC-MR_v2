#!/bin/bash

#SBATCH --chdir=/data/sgg3/liza/SEM_Realv2/scripts # The Working Directory of the job


#Rscript get_tX_UKBB_Automated.R /data/sgg2/liza/SEM_Real/test/data/845.gwas.imputed_v3.both_sexes.tsv.bgz Edu
#Rscript get_tX_UKBB_Automated.R /data/sgg2/liza/SEM_Real/test/data/4080_irnt.gwas.imputed_v3.both_sexes.tsv.bgz SBP
#Rscript get_tX_UKBB_Automated.R /data/sgg2/liza/SEM_Real/test/data/2443.gwas.imputed_v3.both_sexes.tsv.bgz DM
#Rscript get_tX_UKBB_Automated.R /data/sgg2/liza/SEM_Real/test/data/50_irnt.gwas.imputed_v3.both_sexes.tsv.bgz SHeight
#Rscript get_tX_UKBB_Automated.R /data/sgg2/liza/SEM_Real/test/data/20022_irnt.gwas.imputed_v3.both_sexes.tsv.bgz BWeight
#Rscript get_tX_UKBB_Automated.R /data/sgg2/liza/SEM_Real/test/data/2887.gwas.imputed_v3.both_sexes.tsv.bgz PSmoke
Rscript get_tX_UKBB_Automated.R /data/sgg2/liza/SEM_Real/test/data/20002_1111.gwas.imputed_v3.both_sexes.tsv.bgz Asthma
Rscript get_tX_UKBB_Automated.R /data/sgg2/liza/SEM_Real/test/data/20003_1140861958.gwas.imputed_v3.both_sexes.tsv.bgz SVstat
Rscript get_tX_UKBB_Automated.R /data/sgg2/liza/SEM_Real/test/data/20002_1075.gwas.imputed_v3.both_sexes.tsv.bgz MI


Rscript gettingSP_betX.R /data/sgg3/liza/SEM_Realv2/data/ /data/sgg3/liza/SEM_Realv2/data/Asthma_uniq.tsv Asthma
Rscript gettingSP_betX.R /data/sgg3/liza/SEM_Realv2/data/ /data/sgg3/liza/SEM_Realv2/data/SVstat_uniq.tsv SVstat
Rscript gettingSP_betX.R /data/sgg3/liza/SEM_Realv2/data/ /data/sgg3/liza/SEM_Realv2/data/MI_uniq.tsv MI