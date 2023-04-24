################################################################################
# This code paralyzes our reference table during construction for speed
################################################################################
library(igraph)

# Parse in command line output:
cmdArgs <- commandArgs(trailingOnly = TRUE) 
z <- as.character(cmdArgs[1])
partions <- as.numeric(cmdArgs[2])
part <- as.numeric(cmdArgs[3])
d_c <- as.character(cmdArgs[4])


# All summary statistics to be created
request <-  c("steady_num", "steady_str_end_len", "two_concurrency", "num_singletons")

param_all <- c("mu", "rho", "sigma",  "w_0", "w_1", "g_union", "lag") 

source("/home/ot25/Research/JP/HIV/Scripts/string_proc_func.R")
source("/home/ot25/Research/JP/HIV/Scripts/summ_stat_calc.R")


# This function creates the reference table using the cluster
reg_sum_data_function <- function(z,
                                  request,
                                  partions,
                                  part,
                                  d_c){
  
  
  # Setting working directory
  wd_sim <- file.path("/n/scratch3/users/o/ot25/Results/Sum_Discovery",
                      z,
                      "Time_Varying/joint", d_c)
  
  setwd(wd_sim)
  
  # Obtaining all files in the directory 
  all_files_sims <- list.files()
  
  # combining files
  all_files <- c(all_files_sims)
  
  all_files_total <- all_files # This is here so I can recover observed 
  
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
  
  # processing iterations
  iter_base_file <- paste(1,
                          z,
                          "vertexList",
                          base_files, sep = "_")
  
  p_i <- seq(from = 1, 
             to = length(iter_base_file), 
             by = floor(length(iter_base_file)/partions))
  
  if(part < partions){
    iter_base_file <- iter_base_file[p_i[part]:I(p_i[part+1] - 1)]
  }else{
    iter_base_file <- iter_base_file[p_i[part]:num_tp]
  }
  
  # List to store all results
  iter_results <- list()
  
  # Checking if results have already been stored
  wd <- file.path("/home/ot25/Research/JP/HIV/Results/Reference_Table",
                  z,
                  "Raw_Files")
  
  setwd(wd)
  
  file_name <- paste(part, z, "ABC_stat_table.csv", sep = "_")
  ref_files <- list.files()
  
  if(file_name %in% ref_files){
    iter_results[[1]] <- read.csv(file_name)
  }else{
    iter_results[[1]] <- matrix(NA, 0, length(request))
    colnames(iter_results[[1]]) <- request
  }
  
  starter <- dim(iter_results[[1]])[1] + 1
  
  if(starter > length(iter_base_file)){
    print("already done")
    return(NULL)
  }
  
  for(k in starter:length(iter_base_file)){# observed data will go in first row
    # Going back to simulated data
    setwd(wd_sim)
    print(c(k/length(iter_base_file), iter_base_file[k]))
    vertex_search <- grep(paste("^", iter_base_file[k], "_.*$", sep =""),
                          all_files,
                          value = TRUE)
    lag <- str_extract(iter_base_file[k],
                       paste("lag",
                             ".{1}\\d*", sep = "_"))
    
    lag <- as.numeric(str_extract(lag, "-*\\d+"))
    

    G <- list() # stores time dependent graphs
    edge_type_data <- list() # stores edge type information

    for(y in 0:1){ # iterates through time points
      #print(c(k, x, y))
      vertex_file <- grep(paste("tv", y, sep = "_"),
                          vertex_search,
                          value = TRUE)

      vertex_data <- read.csv(vertex_file)

      edge_file <- gsub("vertex",
                        "edge",
                        vertex_file)

      edge_data <- read.csv(edge_file)

      G[[y + 1]] <- make_empty_graph(n = dim(vertex_data)[1],
                                     directed = FALSE)

      V(G[[y + 1]])$name <- vertex_data[, 1]

      edges_1 <- c(t(edge_data))
      edges_1 <- gsub(" ", "", edges_1)
      
      
      G[[y + 1]] <- add_edges(graph = G[[y + 1]],
                              edges = edges_1)
      if(y == 0){ # doesn't exist for second time point
        edge_type_file <- gsub("vertexList",
                               "steady_casual_tracker",
                               vertex_file)
        
        edge_type_data[[y+1]] <- read.csv(edge_type_file, header = TRUE)
        
        # Some single vector matrices are read in wrong
        if(dim(edge_type_data[[y+1]])[2] < 3){
          edge_type_data[[y+1]] <- t(edge_type_data[[y+1]])
        }
        
      }



      # Making graph a simple graph with no self loops or multi-edges
      G[[y + 1]] <- simplify(
        graph = G[[y + 1]],
        remove.multiple = TRUE,
        remove.loops = TRUE,
        edge.attr.comb = igraph_opt("edge.attr.comb") # maps multiedges to one edge
      )

    }

    # Storing results
    temp <- matrix(NA, 1, length(request))
    colnames(temp) <- request

    for(request_i in request){
      # Storing iteration based summary statistics
      temp[1, request_i] <- summ_stat_calc(G = G,
                                           request = request_i,
                                           steady_casual_files = edge_type_data,
                                           lag = lag)
    }
    
    iter_results[[1]] <- rbind(iter_results[[1]], temp)
    
    setwd(wd)
    write.csv(iter_results[[1]], 
              file = file_name, 
              row.names = FALSE)
  }
  
  
  # Setting working directory
  wd <- file.path("/home/ot25/Research/JP/HIV/Results/Reference_Table",
                  z,
                  "Raw_Files")
  
  setwd(wd)
  
  write.csv(iter_results[[1]], 
            file = file_name, 
            row.names = FALSE)
}

reg_sum_data_function(z = z,
                      request = request,
                      partions = partions,
                      part = part,
                      d_c = d_c)


 
