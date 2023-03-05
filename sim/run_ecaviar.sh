source ~/py3/bin/activate
cd ${1}
cd ${2}

python ../../pre_SH.py --dir . --loci ${1} --ite 50

st2=`date +%s.%N`
for i in $(seq 1 50); do ~/utils/CAVIAR-C++/eCAVIAR -o $i -l ../${1}.ld -l ../${1}.ld -z Q$i\.z -z G$i\.z; done
ed2=`date +%s.%N`
echo $(echo $ed2-$st2 | bc -l) > CAVIAR.time

rm *.z
rm *.bse
