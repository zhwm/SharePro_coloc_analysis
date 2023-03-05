module load rstudio-server
cd ${1}
cd ${2}

st2=`date +%s.%N`
Rscript ../../coloc.R
ed2=`date +%s.%N`
echo $(echo $ed2-$st2 | bc -l) > coloc.time

st2=`date +%s.%N`
for i in $(seq 1 50); do Rscript ../../coloc_susie.R ../${1} $i; done
ed2=`date +%s.%N`
echo $(echo $ed2-$st2 | bc -l) > coloc_susie.time

