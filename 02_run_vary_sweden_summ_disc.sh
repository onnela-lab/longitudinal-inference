for d_c in discrete; do

for row in {1..690000..200}; do
			# Submitting the jobs
			echo "${row}, ${d_c}"	
 
			sbatch -o out_Sweden_summ_disc_row_${row}_${d_c}.txt \
			-e e_Sweden_summ_disc_row_${row}_${d_c}.txt \
			--job-name=my_analysis_Sweden_summ_disc_row_${row}_${d_c}.txt \
			02_sweden_summ_disc.sh $row $d_c
			#exit
			sleep 0.6 # pause to be kind to the scheduler	
		done
done
 