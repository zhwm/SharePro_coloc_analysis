#module load StdEnv/2020 r/4.0.0
source ~/py3/bin/activate

cd ${1}
cd ${2}

#for i in 1e-3 1e-4 1e-5 1e-6 1e-7; do Rscript ../../coloc_prior.R $i; done
for i in 1e-3 1e-4 1e-5 1e-6 1e-7; do for j in $(seq 1 50); do python ../../../src/SharePro/sharepro_coloc.py --z Q$j\.fastGWA G$j\.fastGWA --ld ../${1}.ld --sigma $i --save $j\_$i; done; done
