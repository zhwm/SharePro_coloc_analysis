cd ${1}
cd ${2}
st2=`date +%s.%N`
for i in $(seq 1 50); do ~/utils/pwcoco/build/pwcoco --bfile ../${1} --sum_stats1 Q$i\.bse --sum_stats2 G$i\.bse; done
ed2=`date +%s.%N`
echo $(echo $ed2-$st2 | bc -l) > PWCoCo.time

rm *.bse
