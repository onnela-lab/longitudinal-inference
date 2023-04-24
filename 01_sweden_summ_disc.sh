#!/bin/bash
#SBATCH -c 1
#SBATCH -t 0-5 
#SBATCH --mem=3000               
#SBATCH -p short
#SBATCH -e hostname_%j.err   
#SBATCH --mail-type=FAIL 
#SBATCH --mail-user=octavioustalbot@g.harvard.edu

module load gcc/9.2.0 R/4.2.1

# $1 here will be var1
Rscript 01_sweden_sum_discovery.R $1 $2 $3 $4 $5