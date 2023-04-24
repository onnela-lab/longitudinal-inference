partions=200

for d_c in discrete; do
for var1 in Sweden; do
for part in {1..200}; do 

#
echo “${var1}, ${partions}, ${part}, ${d_c}”
#
sbatch -o o_01_reg_sum_data_part_${part}_${d_c}.txt \
-e e_01_reg_sum_data_part_${part}_${d_c}.txt \
--job-name=01_reg_sum_data_part_${part}_${d_c} \
03_reg_sum_data.sh $var1 $partions $part $d_c
#exit
sleep 0.4 # pause to be kind to the scheduler
done
done
done 