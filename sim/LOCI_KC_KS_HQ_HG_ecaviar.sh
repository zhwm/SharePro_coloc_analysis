#!/bin/bash
#
#SBATCH --ntasks 1            # number of tasks
#SBATCH --mem 20G            # memory pool per process
#SBATCH -o LOCI_KC_KS_HQ_HG_ecaviar.out    # STDOUT
#SBATCH -t 06:00:00            # time (D-HH:MM)

./run_ecaviar.sh LOCI KC_KS_HQ_HG
