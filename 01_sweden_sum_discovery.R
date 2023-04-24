###################################################################################################
# This script fixes our set of 5 parameters such that they are comparable to the paper
# or samples free parameters to create the mapping functions. Lag is fixed at 15 for illustration.
# For each row, we vary a single parameter across the parameter space. Then we 
# store the edge list for each iteration so that we can plot at a later date. We fix
# the nodes at 250. 

###################################################################################################
# Parse in command line output:
cmdArgs <- commandArgs(trailingOnly = TRUE) 
#var_change <- as.character(cmdArgs[1])
#var_value <- as.numeric(cmdArgs[2])
#iter_1 <- as.numeric(cmdArgs[3])
row_num <- as.numeric(cmdArgs[1])
iter_1 <- as.numeric(cmdArgs[2])
d_c <- as.character(cmdArgs[3])
var_change <- as.character(cmdArgs[4])
fr_fix <- as.character(cmdArgs[5])

# Reading in our Sweden mechanistic model generator
source("/home/ot25/Research/JP/HIV/Scripts/sweden_mechanisim.R")
setwd("/home/ot25/Research/JP/HIV/Scripts/Sum_Discovery/Sweden")
data <- read.csv("mapping_Sweden_params_1_csv")

g_union <- 1

# Sampling from the prior
sigma <- 1/runif(1, 1, 90)
rho <- 1/runif(1, 1, 50)
w_1 <- 1/runif(1, 1, 61)
w_0 <- 1/runif(1, 1, 40)
mu <- 0.000

if(fr_fix == "fixed"){
  # Changes values if we decide to fix the free parameters
  rho <- 0.3
  sigma <- 0.1
  w_0 <- 0.4
  w_1 <- 0.2
}

# # Iterating through desired parameter
if(var_change == "mu"){
  mu <- 0.000
}else{
  if(var_change == "rho"){
    rho <- data[row_num, "rho"]
  }else{
    if(var_change == "sigma"){
      sigma <- data[row_num, "sigma"]
    }else{
      if(var_change == "w_0"){
        w_0 <- data[row_num, "w_0"]
      }else{
        if(var_change == "w_1"){
          w_1 <- data[row_num, "w_1"]
        }else{
          if(var_change == "g_union"){
            g_union <- g_union
          }
        }
      }
    }
  }
}

# Tracking parameter values
print(c(iter_1, mu, rho, sigma, w_0, w_1))


# Generating data
G <- Sweden_MSM(num_nodes = 250,
                mu = mu,
                rho = rho,
                sigma = sigma,
                w_0 = w_0,
                w_1 = w_1,
                iter = 1000,
                g_union = g_union,
                distribution = d_c,
                lag = 15
)


# Saving all data
## working directory to save files
## tv_0 is time varying previous iteration 
## tv_1 is the time varying final iteration
## E is the edges list of the file
## V is the vertext list of the file

d_c <- paste(d_c, fr_fix, sep = "_")
wd_1 <- file.path("/n/scratch3/users/o/ot25/Results/Sum_Discovery/Sweden/Time_Varying",
                  var_change)
wd_2 <- file.path("/n/scratch3/users/o/ot25/Results/Sum_Discovery/Sweden/Time_Varying",
                  var_change,
                  d_c)

if(!dir.exists(wd_1)){
  dir.create(wd_1)
  dir.create(wd_2)
  setwd(wd_2)
}else{
  if(!dir.exists(wd_2)){
    dir.create(wd_2)
  }
  setwd(wd_2)
}

## Modifying file names
var_value <- data[row_num, var_change]

file_name_0E <- paste(iter_1, "_Sweden_edgeList_", var_change, "_", var_value, "_tv_0.csv",
                      sep = "")
file_name_0V <- paste(iter_1, "_Sweden_vertexList_", var_change, "_", var_value, "_tv_0.csv",
                      sep = "")
file_name_0SC <- paste(iter_1, "_Sweden_steady_casual_tracker_", var_change, "_", var_value, "_tv_0.csv",
                       sep = "")
file_name_1E <- paste(iter_1, "_Sweden_edgeList_", var_change, "_", var_value, "_tv_1.csv",
                      sep = "")
file_name_1V <- paste(iter_1, "_Sweden_vertexList_", var_change, "_", var_value, "_tv_1.csv",
                      sep = "")
file_name_1SC <- paste(iter_1, "_Sweden_steady_casual_tracker_", var_change, "_", var_value, "_tv_1.csv",
                       sep = "")


write.csv(x = get.edgelist(G[[1]]), 
          file = file_name_0E,
          row.names = FALSE)


write.csv(x = matrix(V(G[[1]])$name), 
          file = file_name_0V,
          row.names = FALSE)

write.csv(x = get.edgelist(G[[2]]), 
          file = file_name_1E,
          row.names = FALSE)


write.csv(x = matrix(V(G[[2]])$name), 
          file = file_name_1V,
          row.names = FALSE)

write.csv(x = G[[3]], 
          file = file_name_0SC,
          row.names = FALSE)

print("done")


