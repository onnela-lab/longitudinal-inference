#####################################################################################
# This function creates the parameter matrix from our simulations
#####################################################################################
# Parse in command line output:
cmdArgs <- commandArgs(trailingOnly = TRUE) 
z <- as.character(cmdArgs[1])
d_c <- as.character(cmdArgs[2])

param_all <- c("mu", "rho", "sigma",  "w_0", "w_1", "g_union", "lag") 


source("/home/ot25/Research/JP/HIV/Scripts/string_proc_func.R")
# This function converts the parameter files to parameter data 
# from the population
param_matrix_function <- function(z, d_c){
  
  # z:
  # Population
  # d_c:
  # Discrete vs Continuous

  # Setting working directory
  wd_sim <- file.path("/n/scratch3/users/o/ot25/Results/Sum_Discovery",
                      z,
                      "Time_Varying/joint", d_c)
  
  
  setwd(wd_sim)
  
  # Obtaining all files in the directory 
  all_files <- list.files()
  
  all_files_s_c <- grep("casual", all_files, value = TRUE) 
  
  
  # Limit to only the vertex data
  all_files <- all_files[grep("vertex", all_files)]
  param_value_combo <- list()
  for(j in 1:length(param_all)){
    param_value_combo[[j]] <- joint_processor(str = all_files,
                                              param = param_all[j])
    
    param_value_combo[[j]] <- param_value_combo[[j]][seq(from = 1, to = length(all_files), by = 2)]
  }
  
  # Base files
  base_files <- paste(param_value_combo[[1]],
                      param_value_combo[[2]],
                      param_value_combo[[3]],
                      param_value_combo[[4]],
                      param_value_combo[[5]],
                      param_value_combo[[6]],
                      param_value_combo[[7]], # lag parameter
                      sep = "_")
  
  num_tp <- length(base_files)
  
  param_matrix <- matrix(NA, 
                         num_tp, # first few rows are observed data
                         length(param_all))
  colnames(param_matrix) <- param_all
  
  setwd(wd_sim)
  
  for(j in 1:dim(param_matrix)[1]){  
    for(k in 1:dim(param_matrix)[2]){
      var <- colnames(param_matrix)[k]
      
      temp <- str_extract(base_files[j], # base_files/param_... does not include observed data
                          param_value_combo[[k]][j])
      
      param_matrix[j, k]  <- as.numeric(substr(temp,
                                               nchar(var) + 2,
                                               nchar(temp)))
      
    }
  }
  
  # Setting working directory
  wd <- file.path("/home/ot25/Research/JP/HIV/Results/Reference_Table",
                  z)
  
  setwd(wd)
  
  file_name <- paste(z, "param_table.csv", sep = "_")
  write.csv(param_matrix, file = file_name, row.names = FALSE)
}


param_matrix_function(z, d_c)