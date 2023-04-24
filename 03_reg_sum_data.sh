#!/bin/bash
#SBATCH -c 1
#SBATCH -t 0-12 
#SBATCH --mem=5500                 
#SBATCH -p short
#SBATCH -e hostname_%j.err   
#SBATCH --mail-type=FAIL 
#SBATCH --mail-user=octavioustalbot@g.harvard.edu

module load gcc/9.2.0 R/4.2.1

# $1 here will be var1
Rscript 03_reg_sum_data_function.R $1 $2 $3 $4