
for i in {1..22} X
do
echo chr${i}

time -v sh predict.sh chr${i} 2> resource_usage.chr${i}.log

done
