source ~/py3/bin/activate
cd ${1}
cd ${2}

python ../../pre_SH.py --dir . --loci ${1} --ite 50

st2=`date +%s.%N`
python ~/scratch/SharePro_loc/src/sharepro_loc.py --zld SH.zld --zdir . --N ${3} ${4} --save SH --prefix SH --verbose --K 10
ed2=`date +%s.%N`
echo $(echo $ed2-$st2 | bc -l) > SH6.time

rm *.z
rm Q*.log
rm G*.log
