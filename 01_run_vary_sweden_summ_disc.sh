for d_c in discrete; do
for iter in {1..100}; do
for row_num in {1..30}; do
for var_change in rho sigma w_0 w_1; do
for fr_fix in free fixed; do

echo  ${row_num}, ${iter}, ${d_c}, ${var_change}, ${fr_fix} # printing variable/value/iteration combos
			sbatch -o out_sweden_summ_disc_${row_num}_iter_${iter}_${var_change}_${d_c}_${fr_fix}.txt \
			-e e_sweden_summ_disc_${row_num}_iter_${iter}_${var_change}_${d_c}_${fr_fix}.txt \
			--job-name=my_analysis_sweden_summ_disc_${row_num}_iter_${iter}_${var_change}_${d_c}_${fr_fix} \
			01_sweden_summ_disc.sh $row_num $iter $d_c $var_change $fr_fix
			#exit
			sleep 0.6 # pause to be kind to the scheduler

done
done 
done
done
done