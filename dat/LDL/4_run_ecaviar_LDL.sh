#!/bin/bash
#
#SBATCH --ntasks 1            # number of tasks
#SBATCH --mem-per-cpu 5G            # memory pool per process
#SBATCH -o 4_run_ecaviar_LDL.out    # STDOUT
#SBATCH -t 03:00:00            # time (D-HH:MM)

while read a b c d; do ~/utils/CAVIAR-C++/eCAVIAR -o res_1e-5/ecaviar_LDL_$d -l ss/$d\.ld -l ss/$d\.ld -z ss/LDL.$d\.LDL.txt -z ss/LDL.$d\.$d\.txt; done < LDL.bed
