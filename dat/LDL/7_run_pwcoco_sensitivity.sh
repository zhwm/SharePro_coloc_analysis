#!/bin/bash
#
#SBATCH --ntasks 1            # number of tasks
#SBATCH --mem 10G            # memory pool per process
#SBATCH -o 7_run_pwcoco_sensitivity.out    # STDOUT
#SBATCH -t 3:00:00            # time (D-HH:MM)

for p in 1e-3 1e-4 1e-5 1e-6 1e-7; do while read a b c d; do ~/utils/pwcoco/build/pwcoco --bfile bed/$d --sum_stats1 ss/LDL.$d\.LDL.pwcoco --sum_stats2 ss/LDL.$d\.$d\.pwcoco --coloc_pp 1e-4 1e-4 $p --out res_$p/pwcoco_LDL_$d --maf 0.01 --top_snp 10; done < LDL.bed; done
for p in 1e-3 1e-4 1e-5 1e-6 1e-7; do while read a b c d; do ~/utils/pwcoco/build/pwcoco --bfile bed/$d --sum_stats1 ss/Tg.$d\.Tg.pwcoco --sum_stats2 ss/Tg.$d\.$d\.pwcoco --coloc_pp 1e-4 1e-4 $p --out res_$p/pwcoco_Tg_$d --maf 0.01 --top_snp 10; done < Tg.bed; done
