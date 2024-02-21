#!/bin/bash
#
#SBATCH --ntasks 1            # number of tasks
#SBATCH --mem 10G            # memory pool per process
#SBATCH -o 6_run_csusie_sensitivity.out    # STDOUT
#SBATCH -t 03:00:00            # time (D-HH:MM)

module load r/4.1.2
for p in 1e-3 1e-4 1e-5 1e-6 1e-7; do while read a b c d; do Rscript 6_coloc_susie.R ss/LDL.$d\.LDL.pwcoco ss/LDL.$d\.$d\.pwcoco ss/$d\.ld $p res_$p/csusie_LDL_$d\.txt; done < LDL.bed; done
for p in 1e-3 1e-4 1e-5 1e-6 1e-7; do while read a b c d; do Rscript 6_coloc_susie.R ss/Tg.$d\.Tg.pwcoco ss/Tg.$d\.$d\.pwcoco ss/$d\.ld $p res_$p/csusie_Tg_$d\.txt; done < Tg.bed; done
