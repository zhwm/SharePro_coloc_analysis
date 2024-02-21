#!/bin/bash
#
#SBATCH --ntasks 1            # number of tasks
#SBATCH --mem-per-cpu 5G            # memory pool per process
#SBATCH -o 4_run_ecaviar_Tg.out    # STDOUT
#SBATCH -t 03:00:00            # time (D-HH:MM)

while read a b c d; do ~/utils/CAVIAR-C++/eCAVIAR -o res_1e-5/ecaviar_Tg_$d -l ss/$d\.ld -l ss/$d\.ld -z ss/Tg.$d\.Tg.txt -z ss/Tg.$d\.$d\.txt; done < Tg.bed
