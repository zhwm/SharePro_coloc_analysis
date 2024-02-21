#!/bin/bash
#
#SBATCH --ntasks 40            # number of tasks
#SBATCH --mem-per-cpu 5G            # memory pool per process
#SBATCH -o 3_calculate_ld.out    # STDOUT
#SBATCH -t 03:00:00            # time (D-HH:MM)

module load plink/1.9b_6.21-x86_64
while read a b c d; do plink --bfile bed/$d --extract ss/LDL.$d\.LDL.txt --matrix --r --out ss/$d; done < LDL.bed
while read a b c d; do plink --bfile bed/$d --extract ss/Tg.$d\.Tg.txt --matrix --r --out ss/$d; done < Tg.bed
