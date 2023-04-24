####################################################################################
# This script generates random parameter values from a *reasonable* range 
# to investigate our mechanistic models and the joint distribution of our
# summary statistics
####################################################################################
set.seed(2919865)
setwd("/Volumes/SanDisk/Research/JP/HIV/Data/model_parameters")

rho <- w_1 <- w_0 <- sigma <- round(seq(0.001, 1, length = 30), 3)


Sweden_params <- cbind(rho,
                       sigma,
                       w_0,
                       w_1)

write.csv(Sweden_params,
          "mapping_Sweden_params_1_csv",
          row.names = FALSE)


###############
set.seed(29195)
setwd("/Volumes/SanDisk/Research/JP/HIV/Data/model_parameters")

# The range of possible values is determined using the exponetial distribution and 
# data from the original paper. In later scripts, the numbers are hard coded.

# Sampled parameters
r_quant <- function(x, lamda){
  log(1-x)/(-1*lamda)
}

# Sampling from the prior
m_sigma <- 1/(292.2/30)
q <- r_quant(.9999, m_sigma)
sigma <- rep(1/runif(30000, 1, q),
             23)

m_single <- 1/(163/30)
q <- r_quant(.9999, m_single)
rho <- rep(1/runif(30000, 1, q),
           23)

a_1 = 1/36.3
a_0 = 1/23.1

m_sc <- a_1/sqrt(a_0*0.36+ a_1*(1-0.36))
q <- r_quant(.9999, m_sc)
w_1 <- rep(1/runif(30000, 1, q),
           23)

m_cc <- a_0/sqrt(a_0*0.36+ a_1*(1-0.36)) 
q <- r_quant(.9999, m_cc)
w_0 <- rep(1/runif(30000, 1, q),
           23)

mu <- rep(rep(0.001, 30000),
          23)

g_union <- rep(1, 690000)


lag <- rep(c(-1, 
             seq(5, 70, by = 5),
             seq(80, 150, by = 10)), each = 30000)


Sweden_params <- cbind(mu, 
                       rho,
                       sigma,
                       w_0,
                       w_1,
                       g_union,
                       lag)

Sweden_params <- round(Sweden_params, 10)

write.csv(Sweden_params,
          "matrix_Sweden_params_19_csv",
          row.names = TRUE)

