#!/bin/bash
#
#SBATCH --ntasks 40            # number of tasks
#SBATCH --mem-per-cpu 5G            # memory pool per process
#SBATCH -o 1_get_bed.out    # STDOUT
#SBATCH -t 03:00:00            # time (D-HH:MM)

module load plink/1.9b_6.21-x86_64
while read a b c d; do plink2 --bfile ~/scratch/UKB_Geno/$a --keep ukb_eur_30000.fam --chr $a --from-bp $b --to-bp $c --maf 0.01 --make-bed --out bed/$d; done < gene.bed
