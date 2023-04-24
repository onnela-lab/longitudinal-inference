for var1 in Sweden; do

#
echo “${var1}”
#
sbatch -o o_02_reg_sum_data.txt \
-e e_02_reg_sum_data.txt \
--job-name=02_reg_sum_data \
04_reg_sum_data.sh $var1 
#exit
sleep 0.4 # pause to be kind to the scheduler
done
