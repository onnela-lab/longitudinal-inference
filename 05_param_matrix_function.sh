#!/bin/bash
#SBATCH -c 2
#SBATCH -t 0-3 
#SBATCH --mem=5000                 
#SBATCH -p short
#SBATCH -e hostname_%j.err   
#SBATCH --mail-type=FAIL 
#SBATCH --mail-user=octavioustalbot@g.harvard.edu

module load gcc/9.2.0 R/4.2.1

# $1 here will be var1
Rscript 05_param_matrix_function.R $1 $2