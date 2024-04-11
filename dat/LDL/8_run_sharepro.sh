#!/bin/bash
#
#SBATCH --ntasks 1            # number of tasks
#SBATCH --mem 10G            # memory pool per process
#SBATCH -o 8_run_sharepro.out    # STDOUT
#SBATCH -t 3:00:00            # time (D-HH:MM)

source ~/py3/bin/activate

for p in 1e-3 1e-4 1e-5 1e-6 1e-7; do python ../../src/SharePro/sharepro_coloc.py --z ss/LDL.PCSK9.LDL.pwcoco ss/LDL.PCSK9.PCSK9.pwcoco --ld ss/PCSK9.ld --save res_$p/LDL.PCSK9 --sigma $p; done
for p in 1e-3 1e-4 1e-5 1e-6 1e-7; do python ../../src/SharePro/sharepro_coloc.py --z ss/LDL.APOB.LDL.pwcoco ss/LDL.APOB.APOB.pwcoco --ld ss/APOB.ld --save res_$p/LDL.APOB --sigma $p; done
for p in 1e-3 1e-4 1e-5 1e-6 1e-7; do python ../../src/SharePro/sharepro_coloc.py --z ss/Tg.LPL.Tg.pwcoco ss/Tg.LPL.LPL.pwcoco --ld ss/LPL.ld --save res_$p/Tg.LPL --sigma $p; done
for p in 1e-3 1e-4 1e-5 1e-6 1e-7; do python ../../src/SharePro/sharepro_coloc.py --z ss/Tg.APOC3.Tg.pwcoco ss/Tg.APOC3.APOC3.pwcoco --ld ss/APOC3.ld --save res_$p/Tg.APOC3 --sigma $p; done
for p in 1e-3 1e-4 1e-5 1e-6 1e-7; do python ../../src/SharePro/sharepro_coloc.py --z ss/Tg.ANGPTL3.Tg.pwcoco ss/Tg.ANGPTL3.ANGPTL3.pwcoco --ld ss/ANGPTL3.ld --save res_$p/Tg.ANGPTL3 --sigma $p; done