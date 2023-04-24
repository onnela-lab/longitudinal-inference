#!/bin/bash
#SBATCH -c 2
#SBATCH -t 0-6 
#SBATCH --mem=15000                 
#SBATCH -p short
#SBATCH -e hostname_%j.err   
#SBATCH --mail-type=FAIL 
#SBATCH --mail-user=octavioustalbot@g.harvard.edu

module load gcc/9.2.0 R/4.2.1

# $1 here will be var1
Rscript 04_reg_sum_data_function.R $1