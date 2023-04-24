###################################################################################################
# This script samples our set of 5 parameters such that they are comparable to the paper.
# Next, for each row, we vary a single parameter accross the parameter space. Then we 
# store the edge list for each iteration so that we can plot at a later date. We fix
# the nodes at 250. 

###################################################################################################
# Parse in command line output:
cmdArgs <- commandArgs(trailingOnly = TRUE) 
row_num <- as.numeric(cmdArgs[1])
d_c <- as.character(cmdArgs[2])
 
# Must be installed on the cluster prior 
library(igraph)

# Reading in our Sweden mechanistic model generator
source("/home/ot25/Research/JP/HIV/Scripts/sweden_mechanisim.R")

setwd("/home/ot25/Research/JP/HIV/Scripts/Sum_Discovery/Sweden")
full_data <- read.csv("matrix_Sweden_params_19_csv")

# 10 parameter sets per node
row_nums <- row_num:I(row_num+199)  
subset_data <- full_data[row_nums, ]

for(n_data in 1:dim(subset_data)[1]){ 
  
  # Each row is different
  mu <- 0.000
  rho <- subset_data[n_data, "rho"]
  sigma <- subset_data[n_data, "sigma"]
  w_0 <- subset_data[n_data, "w_0"]
  w_1 <- subset_data[n_data, "w_1"]
  g_union <- 1
  lag <- subset_data[n_data, "lag"]
  iter_1 <- 1
  post_table <- 0
  
  print(paste("complete = ", n_data/dim(subset_data)[1],
              " mu = ", mu, 
              ", rho = ", rho,
              ", sigma = ", sigma,
              ", w_0 = ", w_0,
              ", w_1 = ", w_1,
              ", g_union = ", g_union,
              ", lag = ", lag,
              ", iteration = ", iter_1,
              ", post_table = ", post_table, sep = ""))
  
  G <- Sweden_MSM(num_nodes = 250,
                  mu = mu,
                  rho = rho,
                  sigma = sigma,
                  w_0 = w_0,
                  w_1 = w_1,
                  iter = 1000,
                  g_union = g_union,
                  distribution = d_c,
                  lag = lag
  )
  
  # working directory to save files
  # tv_0 is time varying previous iteration 
  # tv_1 is the time varying final iteration
  # E is the edges list of the file
  # V is the vertex list of the file
  
  if(d_c == "discrete"){
    setwd(paste("/n/scratch3/users/o/ot25/Results/Sum_Discovery/Sweden/Time_Varying/joint/discrete",
                sep = ""))
  }else{
    if(d_c == "continuous"){
      setwd(paste("/n/scratch3/users/o/ot25/Results/Sum_Discovery/Sweden/Time_Varying/joint/continuous",
                  sep = ""))
    }
  }
  
  file_name_0E <- paste(iter_1, "_Sweden_edgeList_mu_",
                        mu,
                        "_rho_", rho,
                        "_sigma_", sigma,
                        "_w_0_", w_0,
                        "_w_1_", w_1,
                        "_g_union_", g_union,
                        "_lag_", lag,
                        "_tv_0.csv",
                        sep = "")
  file_name_0V <- paste(iter_1, "_Sweden_vertexList_mu_",
                        mu,
                        "_rho_", rho,
                        "_sigma_", sigma,
                        "_w_0_", w_0,
                        "_w_1_", w_1,
                        "_g_union_", g_union,
                        "_lag_", lag,
                        "_tv_0.csv",
                        sep = "")
  
  file_name_1E <- paste(iter_1, "_Sweden_edgeList_mu_",
                        mu,
                        "_rho_", rho,
                        "_sigma_", sigma,
                        "_w_0_", w_0,
                        "_w_1_", w_1,
                        "_g_union_", g_union,
                        "_lag_", lag,
                        "_tv_1.csv",
                        sep = "")
  file_name_1V <- paste(iter_1, "_Sweden_vertexList_mu_",
                        mu,
                        "_rho_", rho,
                        "_sigma_", sigma,
                        "_w_0_", w_0,
                        "_w_1_", w_1,
                        "_g_union_", g_union,
                        "_lag_", lag,
                        "_tv_1.csv",
                        sep = "")
  
  file_name_0SC <- paste(iter_1, "_Sweden_steady_casual_tracker_mu_",
                         mu,
                         "_rho_", rho,
                         "_sigma_", sigma,
                         "_w_0_", w_0,
                         "_w_1_", w_1,
                         "_g_union_", g_union,
                         "_lag_", lag,
                         "_tv_0.csv",
                         sep = "")
  
  file_name_1SC <- paste(iter_1, "_Sweden_steady_casual_tracker_mu_",
                         mu,
                         "_rho_", rho,
                         "_sigma_", sigma,
                         "_w_0_", w_0,
                         "_w_1_", w_1,
                         "_g_union_", g_union,
                         "_lag_", lag,
                         "_tv_1.csv",
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
  
  #write.csv(x = G[[4]], 
  #          file = file_name_1SC,
  #          row.names = FALSE)
  
}












