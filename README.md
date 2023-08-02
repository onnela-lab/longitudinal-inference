# longitudinal-inference

This code simulates contacts among men who have sex with men in Stockholm Sweden and demonstrates improved infereces on the mechanism when data are collected in two waves. The improvement is dependent on the timing of data collection and available relationship data demonstrating a need investigate hypothectical population network mechanisms prior to the study. 

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

### Usage

To generate a two network samples of 250 people 15 iterations apart from our discretized network mechanism with mu = 0.01, rho = 0.1, sigma = 0.1, w_0 = .3, w_1 = .2


```
G <- Sweden_MSM(num_nodes = 250,
                mu = 0.01,
                rho = 0.1,
                sigma = 0.1,
                w_0 = .2,
                w_1 = .3,
                iter = 1000,
                g_union = 1,
                distribution = "discrete",
                lag = 15
)
```

We note we are not aggregating edges past the target iteration (g_union = 1).

To calcuate the  summary statistic "propotion os steady relationship maintained between lag points" using a generated network 

```
summary <- summ_stat_calc(G = G,
                          request = "steady_m",
                          steady_casual_files = G[[3]],
                          lag = 15)
```



