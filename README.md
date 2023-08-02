# longitudinal-inference

This code simulates contacts among Men who have sex with Men in Stockholm Sweden and demonstrates improved infereces on the mechanism when data are collected in two waves. The improvement is dependent on the timing of data collection and available relationship data demonstrating a need investigate hypothectical population network mechanisms prior to the study. 

## File Overview

### Core Functions

00_param_gen.R: This script generates random parameter values from a *reasonable* range to investigate our mechanistic models and the joint distribution of our summary statistics

01_sweden_sum_discovery.R: Generates data for our mapping function

02_sweden_sum_discovery.R: Generates our population

03_reg_sum_data_function.R: This code paralyzes our reference table during construction for speed

