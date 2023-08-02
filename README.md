# longitudinal-inference

This code simulates contacts among Men who have sex with Men in Stockholm Sweden and demonstrates improved infereces on the mechanism when data are collected in two waves. The improvement is dependent on the timing of data collection and available relationship data demonstrating a need investigate hypothectical population network mechanisms prior to the study. 

## File Overview

### Data Generation

00_param_gen.R: This script generates random parameter values from a *reasonable* range to investigate our mechanistic models and the joint distribution of our summary statistics.

01_sweden_sum_discovery.R: Generates data for our mapping function.

02_sweden_sum_discovery.R: Generates our population.

### Data Aggregation/Paralyzation

03_reg_sum_data_function.R: This code paralyzes our reference table for Approximate Bayesian Computation during construction for speed.

04_reg_sum_data_function.R: This script combines all raw files into the reference table.

05_param_matrix_function.R: This function recreates the parameter matrix from our simulations where the parameters are in the file name.

### Figures
06_summ_discovery_plots.R: This script calculates our summary statistic plots.

07_figures.R: This script generates the lag plots for our paper





