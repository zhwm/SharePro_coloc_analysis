module load rstudio-server
cd ${1}
mkdir ${2}\_${3}\_${4}\_${5}
cd ${2}\_${3}\_${4}\_${5}
Rscript ../../simulate_traits.R ../${1} ${2} ${3} ${4} ${5}
for i in $(seq 1 50); do ~/utils/gcta_1.93.2beta/gcta64 --bfile ../${1} --pheno Q.phen --mpheno $i --fastGWA-lr --out Q$i; done
for i in $(seq 1 50); do ~/utils/gcta_1.93.2beta/gcta64 --bfile ../${1} --pheno G.phen --mpheno $i --fastGWA-lr --out G$i; done
