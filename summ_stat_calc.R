################################################################################################
# This function calculates different summary statistics depending on the input.
# t_varying should be TRUE for any graph inputed in this algorithm.

################################################################################################
source("/home/ot25/Research/JP/HIV/Scripts/common_edges.R")

summ_stat_calc <- function(G, 
                           request,
                           steady_casual_files = NULL,
                           lag){
  # G:
  # Graph we are obtaining the summary statistic from
  # request:
  # Summary statistic we are requesting
  # num_singletons: number of nodes with degree 0
  # two_concurrency: mean percentage of people in steady relationships also in casual ones
  # steady_str_end_len: average steady relationship length of relationships that start and end during the study
  # steady_num: proportion of relationships that are steady
  
  # steady_casual_file:
  # For the sweden model, these are the accompanying files that contains information on whether 
  # the edges is casual or steady
  
  
  if(length(G) == 2){
    g_0 <- G[[1]]
    g_1 <- G[[2]]
    
    E_0 <- get.edgelist(g_0)
    E_1 <- get.edgelist(g_1)
  }else{
    g_1 <- G
    
    E_1 <- get.edgelist(g_1)
  }
  
  last_iteration <- max(steady_casual_files[[1]][, "iteration"])
  
  # lag -1 implies same time point. 
  lag <- lag + 1
  
  if(request == "num_singletons"){
    sum_stat_0 <- sum(degree(g_0) == 0)
    sum_stat_1 <- sum(degree(g_1) == 0)
    
    sum_stat_0 <- sum_stat_0/length(V(g_0))
    sum_stat_1 <- sum_stat_1/length(V(g_1))
    
    sum_stat <- (sum_stat_0 + sum_stat_1)/2
  }else{
    if(request == "steady_str_end_len"){
      if(is.null(steady_casual_files)){
        sum_stat <- NA 
        
      }else{
        initial_steady <- steady_casual_files[[1]][which(steady_casual_files[[1]][, "type"]
                                                         == "steady" | 
                                                           steady_casual_files[[1]][, "type"]
                                                         == "s"), ,drop = FALSE]
        
        # which relationships were ongoing during duration of study and ends
        # during duration of the study
        nn <- which(initial_steady[, "start"] >= last_iteration - lag  - 12 & # added 3 iterations
                      initial_steady[, "end"] <= last_iteration)
        
        if(length(nn) > 0){
          duration <- initial_steady[nn, "end"] - initial_steady[nn, "start"]
          sum_stat <- mean(duration, na.rm = TRUE)
        }else{
          sum_stat <- NA
        }  
      }
    }else{
      if(request == "two_concurrency"){
        if(is.null(steady_casual_files)){
          sum_stat <- NA
          
        }
        end_steady_0 <- steady_casual_files[[1]][which(steady_casual_files[[1]][, "type"]
                                                       == "steady" | 
                                                         steady_casual_files[[1]][, "type"]
                                                       == "s"), ,drop = FALSE]
        
        nn <- which(end_steady_0[, "start"] <= last_iteration - lag & # steady edges present at first lag
                      end_steady_0[, "end"] > last_iteration - lag) 
        
        end_steady_0 <- end_steady_0[nn, ]
        st_0 <- c(end_steady_0[, "e1"], end_steady_0[, "e2"])
        
        end_casual_0 <- steady_casual_files[[1]][which(steady_casual_files[[1]][, "type"]
                                                       == "casual" | 
                                                         steady_casual_files[[1]][, "type"]
                                                       == "c"), ,drop = FALSE]
        
        nn <- which(end_casual_0[, "start"] == last_iteration - lag) # casual edges are deleted day they are added
        end_casual_0 <- end_casual_0[nn, ]
        
        cs_0 <- c(end_casual_0[, "e1"], end_casual_0[, "e2"])
        
        if(length(st_0) != 0){
          sum_stat_0 <- mean(cs_0 %in% st_0)
        }else{
          sum_stat_0 = 0
        }
        
        
        end_steady_1 <- steady_casual_files[[1]][which(steady_casual_files[[1]][, "type"]
                                                       == "steady" | 
                                                         steady_casual_files[[1]][, "type"]
                                                       == "s"), ,drop = FALSE]
        
        nn <- which(end_steady_1[, "end"] > last_iteration) # edges present at last iteration
        
        end_steady_1 <- end_steady_1[nn, ]
        
        st_1 <- c(end_steady_1[, "e1"], end_steady_1[, "e2"])
        
        end_casual_1 <- steady_casual_files[[1]][which(steady_casual_files[[1]][, "type"]
                                                       == "casual" | 
                                                         steady_casual_files[[1]][, "type"]
                                                       == "c"), ,drop = FALSE]
        
        nn <- which(end_casual_0[, "start"] == last_iteration) # casual edges are deleted day they are added
        end_casual_0 <- end_casual_0[nn, ]
        
        cs_1 <- c(end_casual_1[, "e1"], end_casual_1[, "e2"])
        
        if(length(st_1) != 0){
          sum_stat_1 <- mean(cs_1 %in% st_1)
        }else{
          sum_stat_1 = 0
        }
        
        sum_stat <- (sum_stat_0 + sum_stat_1)/2
        
      }else{
        if(request == "steady_num"){
          end_steady_0 <- steady_casual_files[[1]][which(steady_casual_files[[1]][, "type"]
                                                         == "steady" | 
                                                           steady_casual_files[[1]][, "type"]
                                                         == "s"), ,drop = FALSE]
          
          nn_0 <- which(end_steady_0[, "start"] <= last_iteration - lag & 
                          end_steady_0[, "end"] > last_iteration - lag)
          
          sum_stat_0 <- length(nn_0)
          
          nn_1 <- which(end_steady_0[, "end"] > last_iteration)
          
          sum_stat_1 <- length(nn_1)
          
          sum_stat_0 <- sum_stat_0/dim(E_0)[1]
          sum_stat_1 <- sum_stat_1/dim(E_1)[1]
          
          
          sum_stat <- (sum_stat_0 + sum_stat_1)/2
        }
      }
    }
  }
  
  
  return(sum_stat) 
}
