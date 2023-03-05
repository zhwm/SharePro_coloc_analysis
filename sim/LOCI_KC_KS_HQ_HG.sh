#!/bin/bash
#
#SBATCH --ntasks 1            # number of tasks
#SBATCH --mem 20G            # memory pool per process
#SBATCH -o LOCI_KC_KS_HQ_HG.out    # STDOUT
#SBATCH -t 06:00:00            # time (D-HH:MM)

./simulate_traits.sh LOCI KC KS HQ HG
./run_SH.sh LOCI KC_KS_HQ_HG 1000 25000
./run_PWCoCo.sh LOCI KC_KS_HQ_HG
./run_coloc.sh LOCI KC_KS_HQ_HG
