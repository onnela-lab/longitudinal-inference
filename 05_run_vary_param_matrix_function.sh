for d_c in discrete; do
for var1 in Sweden; do

#
echo “${var1}, ${d_c}”
#
sbatch -o o_param_matrix_${d_c}.txt \
-e e_param_matrix_${d_c}.txt \
--job-name=param_matrix_${d_c} \
05_param_matrix_function.sh $var1 $d_c
#exit
sleep 0.4 # pause to be kind to the scheduler
done
done