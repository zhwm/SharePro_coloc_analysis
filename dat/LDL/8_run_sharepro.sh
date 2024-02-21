#!/bin/bash
#
#SBATCH --ntasks 1            # number of tasks
#SBATCH --mem 10G            # memory pool per process
#SBATCH -o 8_run_sharepro.out    # STDOUT
#SBATCH -t 3:00:00            # time (D-HH:MM)

source ~/py3/bin/activate

for p in 1e-3 1e-4 1e-5 1e-6 1e-7; do python ../../src/sharepro_bse.py --z ss/LDL.PCSK9.LDL.pwcoco ss/LDL.PCSK9.PCSK9.pwcoco --ld ss/PCSK9.ld ss/PCSK9.ld --save res_$p/SH_LDL.PCSK9.txt --K 10 --sigma $p --verbose; done
for p in 1e-3 1e-4 1e-5 1e-6 1e-7; do python ../../src/sharepro_bse.py --z ss/LDL.APOB.LDL.pwcoco ss/LDL.APOB.APOB.pwcoco --ld ss/APOB.ld ss/APOB.ld --save res_$p/SH_LDL.APOB.txt --K 10 --sigma $p --verbose; done
for p in 1e-3 1e-4 1e-5 1e-6 1e-7; do python ../../src/sharepro_bse.py --z ss/Tg.LPL.Tg.pwcoco ss/Tg.LPL.LPL.pwcoco --ld ss/LPL.ld ss/LPL.ld --save res_$p/SH_Tg.LPL.txt --K 10 --sigma $p --verbose; done
for p in 1e-3 1e-4 1e-5 1e-6 1e-7; do python ../../src/sharepro_bse.py --z ss/Tg.APOC3.Tg.pwcoco ss/Tg.APOC3.APOC3.pwcoco --ld ss/APOC3.ld ss/APOC3.ld --save res_$p/SH_Tg.APOC3.txt --K 10 --sigma $p --verbose; done
for p in 1e-3 1e-4 1e-5 1e-6 1e-7; do python ../../src/sharepro_bse.py --z ss/Tg.ANGPTL3.Tg.pwcoco ss/Tg.ANGPTL3.ANGPTL3.pwcoco --ld ss/ANGPTL3.ld ss/ANGPTL3.ld --save res_$p/SH_Tg.ANGPTL3.txt --K 10 --sigma $p --verbose; done
