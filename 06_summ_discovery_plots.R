################################################################################
# This script calculates our summary statistic plots up top.
# The joint distribution of our parameters
# and summary statistics are also available towards the bottom.

################################################################################
library(reshape)
library(stringr)
library(igraph)
library(ggplot2)
library(matlib)
library(glmnet)
library(gridExtra)
library(ica)
library(RColorBrewer)
library(grDevices)
library(tidyverse)

# Set parameters you are creating plots for (ordering matters)
param_all <- c("mu", "rho", "sigma",  "w_0", "w_1", "g_union", "lag") # Sweden parameters

# Iterations to check
g_union_check = c(1)

# This function calculates our summary statistic of choice
source("/Volumes/SanDisk/Research/JP/HIV/Scripts/summ_stat_calc.R")
source("/Volumes/SanDisk/Research/JP/HIV/Scripts/sweden_mechanisim.R")

# Population to consider (keep in same order as working directory)
population <- c("Sweden")

request <- c("num_singletons", "steady_str_end_len",
             "two_concurrency", "steady_num")

dist_types <- "No_Selection"

quant <- c(0.01)

# This function creates the parameter specific heat maps
multi_heat <- function(observed_data, 
                       simulated_data,
                       euclidean,
                       param_all,
                       stat_table,
                       dist_types,
                       g_union_check,
                       lag_m,
                       z,
                       bins,
                       plots = TRUE,
                       prior = FALSE,
                       reg_adj, 
                       root = TRUE){
  
  combos <- combn(param_all[1:5], 2) # getting all pairwise combos
  
  # initializing observed parameter discrepancy data
  if(euclidean){
    MRSSE <- matrix(NA, 1, 1)
    MRSSE <- MRSSE[-1, 1, drop = FALSE]
    colnames(MRSSE) <- paste(param_all[1:5], collapse = "_")
  }else{
    MRSSE <- matrix(NA, 1, 6)
    MRSSE <- MRSSE[-1, , drop = FALSE]
    
    colnames(MRSSE) <- c(paste("_", param_all[1:5], sep = ""), "all")
  }
  
  if(stat_table){ # in case data was generated previously
    setwd("/Volumes/SanDisk/Research/JP/HIV/Data/ref_table")
    file_name <- paste(z, "Euclidean", euclidean ,paste("ABC_stat_MRSSE_", dist_types,
                                                        euclidean ,"_table.csv", sep = ""), sep = "_")
    final_d <- read.csv(file_name)
  }else{
    setwd("/Volumes/SanDisk/Research/JP/HIV/Data/ref_table")
    
    for(i in 1:dim(observed_data)[1]){
      print(i/dim(observed_data)[1])
      # Combining observed data with sample of simulated data
      full_data_0 <- rbind(observed_data[i, ], 
                           simulated_data)
      
      full_data_1 <- final_function(param_matrix = full_data_0[, param_all, drop = FALSE],
                                    reg_sum_data = full_data_0[, -which(colnames(full_data_0) %in% param_all), drop = FALSE],
                                    dist_types = dist_types,
                                    g_union_check = g_union_check,
                                    lag = lag_m,
                                    z = z)
      
      
      full_data <- full_data_1[[1]]
      
      
      # Identifying posterior data
      full_sim_data <- full_data
      full_dist_data <- full_sim_data[-1, dist_types]
      quant_threshold <- quantile(full_dist_data, probs = quant, na.rm = TRUE)
      
      mean_p <- NULL
      sd_p <- NULL
      for(i_p in 1:length(full_sim_data[1:5])){ #Normalizing parameter data for the discrepancy
        mean_p[i_p] <- mean(full_sim_data[, param_all[i_p]])
        sd_p[i_p] <- sd(full_sim_data[, param_all[i_p]])
        
        if("mu" %in% param_all & # when mu is fixed to be constant dont regress on it
           i_p == which(param_all == "mu")){
          full_sim_data[, param_all[i_p]] <- 0.000
        }else{
          if(sd_p[i_p] == 0){
            full_sim_data[, param_all[i_p]] <- mean_p[i_p]
          }else{
            full_sim_data[, param_all[i_p]] <- (full_sim_data[, param_all[i_p]] - mean_p[i_p])/sd_p[i_p]
          }
        }
      }
      
      # first line is observed data
      if(prior){
        post_data <- full_sim_data
      }else{
        post_data <- full_sim_data[c(1, which(full_dist_data <= quant_threshold) + 1), 
                                   ,drop = FALSE ] #+1 resets order
      }
      
      if(reg_adj){
        # removing observed data
        observed_data_post <- post_data[1, ,drop = FALSE]
        post_data <- post_data[-1, , drop = FALSE]
        
        temp_request <- colnames(simulated_data)[-which(colnames(simulated_data) %in% param_all)]
        # Adjusting posterior
        for(k in 1:length(param_all)){
          #print(param_all[k])
          formula <- paste(param_all[k], "~", paste(temp_request, collapse = "+"))
          lm1 <- lm(formula = formula, data = as.data.frame(post_data))
          
          # observed estimate
          obs_est <- predict(lm1, data.frame(full_sim_data[1, temp_request, drop = FALSE]))
          
          # posterior residuals (neighborhood within s_obs)
          post_resid <- resid(lm1)
          
          # Adjustment
          theata_c <- obs_est + post_resid
          
          # Making negative predicted values 0
          theata_c <- ifelse(theata_c < 0, 0, theata_c)
          
          # sub
          post_data[, param_all[k]] <- theata_c
        }
        
        post_data <- rbind(observed_data_post,
                           post_data)
        
      }
  
      
      # Calculating RSSE
      if(euclidean){
        RSSE <- NULL
      }else{
        RSSE <- matrix(NA, 0, 6)
        colnames(RSSE) <- c(paste("_", param_all[1:5], sep = ""), "all")
      }
      
      for(j in 2:dim(post_data)[1]){ # first line is observed data
        
        if(euclidean){
          RSSE[j - 1] <- norm(post_data[j, param_all[1:5]] - post_data[1, param_all[1:5]],
                              type = "2")^2
          
          if(j == dim(post_data)[1]){
            RSSE <- as.matrix(RSSE)
            colnames(RSSE) <- colnames(MRSSE)
          }
        }else{
          # first line is observed data
          
          RSSE_1 <- matrix(NA, 1, 6)
          for(k in 1:dim(RSSE)[2]){
            if(k == 6){
              RSSE_1[1, k] <- norm(post_data[j, param_all[1:5]] - post_data[1, param_all[1:5]],
                                   type = "2")^2
            }else{
              RSSE_1[1, k] <- norm(post_data[j, param_all[k]] - 
                                     post_data[1, param_all[k]],
                                   type = "2")^2
            }
          }
          
          colnames(RSSE_1) <- c(param_all[1:5], "all")
          
          RSSE <- rbind(RSSE_1, RSSE)
          
        }
        
      }
      
      if(root){
        MRSSE_1 <- t(as.matrix(sqrt(apply(RSSE, MARGIN = 2, FUN = mean)))) 
      }else{
        MRSSE_1 <- t(as.matrix(apply(RSSE, MARGIN = 2, FUN = mean))) 
      }

      MRSSE <- rbind(MRSSE, MRSSE_1)
      
    }
    
    
    
    final_d <- cbind(observed_data[, param_all[1:5]], MRSSE)
    
    file_name <- paste(z, "Euclidean", euclidean ,paste("ABC_stat_MRSSE_", dist_types,
                                                        euclidean ,"_table.csv", sep = ""), sep = "_")
    write.csv(final_d, 
              file_name,
              row.names = FALSE)
    
  }
  
  return(final_d[, colnames(MRSSE)]) # return average discrepancy
}
# This function helps eliminate the iteration and the 
# time varying notation from the string names
# the goal is to automatically determine how many numbers
# we are iterating though for our parameter
character_remove <- function(all_strings){
  # all_stings:
  # strings you seek to remove. Must be of paticular format
  
  all_strings_1 <- all_strings
  
  # stage 1 iter remove
  iter_removed_1 <- unique(substr(all_strings, 3, nchar(all_strings) - 9))
  iter_removed_1 <- ifelse(substring(iter_removed_1,1,1) == 0, 
                           substring(iter_removed_1,2,nchar(iter_removed_1)),
                           iter_removed_1)
  iter_removed_1 <- ifelse(substring(iter_removed_1,1,1) == "_", 
                           substring(iter_removed_1,2,nchar(iter_removed_1)),
                           iter_removed_1)
  
  # If we go to 3 digit iterations this will no longer work
  # for(j in 1:length(iter_removed_1)){ # removing potention "_" in beginning of string
  #   if(substr(iter_removed_1[j], 1, 1) == "_"){
  #     iter_removed_1[j] <- substr(iter_removed_1[j],
  #                                 2,
  #                                 nchar(iter_removed_1[j]))
  #   }
  # }
  
  # Just seraching for edgelist lists
  iter_removed_edge <- unique(grep("edgeList", iter_removed_1, value = TRUE))
  
  # returning 1 iteration of all of the parameters independent of time varying
  return(iter_removed_edge)
  
}

# This function obtains the value paramter from the string
processor <- function(all_strings){
  all_strings_1 <- all_strings
  
  all_strings_2 <- str_replace(all_strings_1, 
                               "^.*([0-9]+\\.+[0-9]*+).*$", 
                               "\\1")
  
  # If the all_strings_1 wasn't modified, its because there were
  # no numbers with decimal points
  if(all(all_strings_1 == all_strings_2)){
    all_strings_2 <- str_replace(all_strings_2, 
                                 "^.*_([0-9]+)$", 
                                 "\\1")
  }
  
  # Identifying numbers that were not converted due to no
  # "." being present between numbers
  n_0 <- suppressWarnings(which(is.na(as.numeric(all_strings_2))))
  
  if(length(n_0) > 0){
    all_strings_2[n_0] <- str_replace(all_strings_2[n_0], 
                                      "^.*([0-9]+).*$", 
                                      "\\1")
  }
  
  # Getting the index order to retrun
  n_order <- order(as.numeric(all_strings_2))
  all_strings_2 <- sort(as.numeric(all_strings_2))
  return(list(all_strings_2, n_order))
}


# This function creates distances and outputs 
# final matrix
final_function <- function(param_matrix,
                           reg_sum_data,
                           dist_types,
                           g_union_check,
                           lag, 
                           z){
  
  nn <- which(complete.cases(reg_sum_data))
  reg_sum_data <- reg_sum_data[nn, ]
  param_matrix <- param_matrix[nn, ]
  
  
  
  # Standardized summary data
  std_sum_data <- scale(reg_sum_data)
  
  
  # Creating the distance columns (different metrics will be used)
  dist_matrix <- matrix(NA, dim(std_sum_data)[1], length(dist_types))
  colnames(dist_matrix) <- dist_types
  
  # Contribution looks at which variables contributes to distance 
  # (percentage explained as well in projections)
  
  L2_calc <- function(x){
    # calculating L2 norm
    
    # x
    # row of the summary matrix
    # 1st row is observed data
    x <- x - std_sum_data[1, ]
    L2 <- norm(x, type = "2")
    
    return(L2)
    
  }
  
  dist_matrix[, dist_types] <- apply(std_sum_data, # use standardized data for L2
                                        MARGIN = 1, 
                                        FUN = L2_calc) 
  
  contribution <- request
  
  final <- cbind(param_matrix, 
                 reg_sum_data)
  final <- cbind(final, dist_matrix)
  
  return(list(final, 
              request))
  
}

# This function creates our ABC posterior
ABC_ref <- function(request,
                    population,
                    param_all,
                    num_samples,
                    dist_types,
                    quant,
                    g_union_check,
                    stat_table = FALSE,
                    reg_adj = FALSE,
                    return_grid = FALSE,
                    cheat = TRUE,
                    three = FALSE,
                    lag = 0,
                    prior = FALSE){
  # request:
    # Summary statistics we are requesting
      # num_edges: number of edges
      # edges_m: proportion of edges maintained between time points
      # num_tri: number of triangles in the graph
      # quant_25: 25th degree distribution percentile
      # quant_50: 50th degree distribution percentile
      # quant_75: 75th degree distribution percentile
      # quant_max: maximum of the degree distribution 
      # quant_mean: mean of the degree distribution 
      # median_local_cc: median of the local clustering coefficient
      # mean_local_cc: mean of the local clustering coefficient 
      # num_con_comp: number of components with 2 or more nodes
      # size_LCC: size of the largest connected component
      # num_singletons: number of nodes with degree 0
      # num_exclusive_pairs: number of monogamous pairs
      # prolonged_single:
      # entering
  # population:
    # Country/City we are looking at
      # Sweden
      # Amsterdam
  # param_all
    # parameters to be considered (each model will be different)
  # num_samples
    # number of samples we are considering from the prior distribution
  # dist_types:
    # Distances to be included in the ABC reference table
      # L2
        # L2 norm
  # quant:
    # percentile to threshold the distance metric
  # g_union_check:
    # which union iteration to check
  # stat_table:
    # TRUE: we already calcualted the stat table
    # FALSE: we still need to calculate stat table
  # reg_adj:
    # Regression adjustment
  # return_grid
    # Return ggplots from posterior
  # three:
    # If true calculate 3 different observed versions of the graph
    # results
  # prior:
    # Wheter you investigate the prior instead of the posterior
  
  for(z in population){
    
    setwd("/Volumes/SanDisk/Research/JP/HIV/Data/ref_table")
    # Generating parameter matrix
    file <- paste(z, "param_table.csv", sep = "_")
    param_matrix <- read.csv(file)
    
    
    # Generating reference table
    # Saving time and reading long portion of data in
    
    setwd("/Volumes/SanDisk/Research/JP/HIV/Data/ref_table")
    file_name <- paste(z, "ABC_stat_table.csv", sep = "_")
    
    reg_sum_data <- read.csv(file_name)
    reg_sum_data <- reg_sum_data[, request, drop = FALSE] # removing extra column

    # observed data
    o_data <- cbind(param_matrix[1,,drop = FALSE], reg_sum_data[1,, drop = FALSE])
    
    for(i_o in 1:dim(o_data)[1]){
      param_matrix[1, param_all] <- as.matrix(o_data[i_o, param_all])
      reg_sum_data[1, request] <- as.matrix(o_data[i_o, request])
      
      
      
      final_full <- final_function(param_matrix = param_matrix,
                                   reg_sum_data = reg_sum_data,
                                   dist_types = dist_types,
                                   g_union_check = g_union_check,
                                   lag = lag,
                                   z = z)
      
      final <- final_full[[1]]
      contribution <- final_full[[2]]
      
      # Creating Posterior Plots
      post_plots <- function(full_data, dist_k, quant_x, prior = FALSE){
        # full_data:
        # data to be used to obtain the posterior
        # dist_k:
        # particular distances to use when creating the posterior 
        # plot
        # quant_x:  
        # percentile to check when selecting samples
        
        
        # identifying data
        full_sim_data <- full_data
        full_dist_data <- full_sim_data[-1, dist_k]
        quant_threshold <- quantile(full_dist_data, probs = quant_x, na.rm = TRUE)
        print(quant_threshold)
        if(prior){
          post_data <- full_sim_data
        }else{
          post_data <- full_sim_data[which(full_dist_data <= quant_threshold) + 1, ,drop = FALSE ] #+1 resets order
        }
        
        # storing the plots
        p_length <- 5
        
        all_ggplots = NULL
        for(k in 1:p_length){ # iterate through the 5 parameters
          if(param_all[k] == "gamma"){
            x_lim = c(1.5, 4)
            xlab = expression(gamma)
          }else{
            if(param_all[k] == "p_0"){
              x_lim = c(0, 0.8)
              xlab = param_all[k] 
            }else{
              if(param_all[k] == "p_cas"){
                x_lim = c(0, 1)
                xlab = param_all[k]
              }else{
                if(param_all[k] == "k_max"){
                  x_lim = c(5, 100)
                  xlab = param_all[k]
                }else{
                  if(param_all[k] == "migration"){
                    x_lim = c(0, 0.5)
                    xlab = param_all[k]
                  }else{
                    if(param_all[k] == "mu"){
                      x_lim = c(0, 0.5)
                      xlab = expression(mu)
                    }else{
                      if(param_all[k] == "rho"){
                        x_lim = c(0, 1)
                        xlab = expression(rho)
                      }else{
                        if(param_all[k] == "sigma"){
                          x_lim = c(1/60, 1/2)
                          xlab = expression(sigma)
                        }else{
                          if(param_all[k] == "w_0" & z == "Sweden"){
                            x_lim = c(0, 1)
                            xlab = param_all[k]
                          }else{
                            if(param_all[k] == "w_1"){
                              x_lim = c(0, 1)
                              xlab = param_all[k]
                            }else{
                              if(param_all[k] == "g_union"){
                                x_lim <- c(1, 12)
                              }else{
                                if(param_all[k] == "theta"){
                                  x_lim <- c(0, 1)
                                  xlab = param_all[k]
                                }else{
                                  if(param_all[k] == "m"){
                                    x_lim <- c(10, 60)
                                    xlab = "mean"
                                  }else{
                                    if(param_all[k] %in% c("p_r", "p_d")){
                                      x_lim <- c(0, 1)
                                      xlab <- param_all[k]
                                    }else{
                                      if(param_all[k] %in% c("w_0", "delta") & z == "PRL"){
                                        x_lim <- c(0, 2)
                                        xlab <- param_all[k]
                                      }else{
                                        if(param_all[k] %in% c("size")){
                                          x_lim <- c(140, 160)
                                          xlab <- param_all[k]
                                        }else{
                                          if(param_all[k] %in% c("theta")){
                                            x_lim <- c(0.1, 0.3)
                                            xlab <- param_all[k]
                                          }else{
                                            if(param_all[k] %in% c("std")){
                                              x_lim <- c(2, 4)
                                              xlab <- param_all[k]
                                            }
                                          }
                                        }
                                      }
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
          
          
          if(k == 1){
            y_lab = "Count"
          }else{
            y_lab = ""
          }
          
          # Making the plots
          
          all_ggplots[[k]] <- p <- ggplot(as.data.frame(post_data), aes_string(x = param_all[k]))+
            theme_bw() +
            geom_histogram(aes(y = stat(count)), bins = 20, 
                           color = "black", fill = "grey") +
            #scale_y_continuous(labels = scales::percent) + 
            
            xlab(xlab)+
            ylab(y_lab) +
            
            scale_x_continuous(limits = x_lim, oob = scales::oob_keep) + 
            ylim(c(0, dim(post_data)[1])) + 
            geom_vline(xintercept = full_data[1, param_all[k]], col = "red", linetype="dotted") + #observed
            geom_vline(xintercept = quantile(post_data[, param_all[k]], 
                                             probs = 0.05, 
                                             na.rm = TRUE),
                       col="blue") + 
            geom_vline(xintercept = quantile(post_data[, param_all[k]], 
                                             probs = 0.95, 
                                             na.rm = TRUE),
                       col="blue") 
          
          #return(all_plots)
        }
        
        
        if(return_grid &! reg_adj){
          return(all_ggplots)
        }else{
          if(length(param_all) != 1 & !reg_adj){
            grid.arrange(all_ggplots[[1]],
                         all_ggplots[[2]],
                         all_ggplots[[3]],
                         all_ggplots[[4]],
                         all_ggplots[[5]],
                         ncol = 5,
                         top = paste(z, " Posterior", "--", 
                                     quant_x*100, "th percentile \n", dist_k, "\n (",
                                     g_union_check, " Iteration(s), ",
                                     "lag = ", lag, ")", sep = ""))
            
            return(all_ggplots)
          }else{
            
            if(!reg_adj){
              grid.arrange(all_ggplots[[1]],
                           ncol = 1,
                           top = paste(z, " Posterior", "--", 
                                       quant_x*100, "th percentile \n", dist_k, "\n (",
                                       g_union_check, " Iteration(s), ",
                                       "lag = ", lag, ")", sep = ""))
            }
          }
          
        }
        
        
        
        
        
        
        # Regression adjustment
        if(reg_adj){
          
          # Adjusting posterior
          for(k in 1:length(param_all)){
            #print(param_all[k])
            formula <- paste(param_all[k], "~", paste(request, collapse = "+"))
            lm1 <- lm(formula = formula, data = as.data.frame(post_data))
            
            # observed estimate
            obs_est <- predict(lm1, data.frame(full_data[1, request, drop = FALSE]))
            
            # posterior residuals (neighborhood within s_obs)
            post_resid <- resid(lm1)
            
            # Adjustment
            theata_c <- obs_est + post_resid
            
            # Making negative predicted values 0
            theata_c <- ifelse(theata_c < 0, 0, theata_c)
            
            # sub
            post_data[, param_all[k]] <- theata_c
          }
          
          # Making plots
          
          # storing the plots
          all_ggplots = NULL
          for(k in 1:length(param_all[1:5])){
            if(param_all[k] == "gamma"){
              x_lim = c(1.5, 4)
              xlab = expression(gamma)
            }else{
              if(param_all[k] == "p_0"){
                x_lim = c(0, 0.8)
                xlab = param_all[k] 
              }else{
                if(param_all[k] == "p_cas"){
                  x_lim = c(0, 1)
                  xlab = param_all[k]
                }else{
                  if(param_all[k] == "k_max"){
                    x_lim = c(5, 100)
                    xlab = param_all[k]
                  }else{
                    if(param_all[k] == "migration"){
                      x_lim = c(0, 0.2)
                      xlab = param_all[k]
                    }else{
                      if(param_all[k] == "mu"){
                        x_lim = c(0, 0.5)
                        xlab = expression(mu)
                      }else{
                        if(param_all[k] == "rho"){
                          x_lim = c(0, 1)
                          xlab = expression(rho)
                        }else{
                          if(param_all[k] == "sigma"){
                            x_lim = c(0, 0.5)
                            xlab = expression(sigma)
                          }else{
                            if(param_all[k] == "w_0"){
                              x_lim = c(0, 1)
                              xlab = param_all[k]
                            }else{
                              if(param_all[k] == "w_1"){
                                x_lim = c(0, 1)
                                xlab = param_all[k]
                              }else{
                                if(param_all[k] == "theta"){
                                  x_lim <- c(0, 1)
                                  xlab = param_all[k]
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
            
            
            if(k == 1){
              y_lab = "Count"
            }else{
              y_lab = ""
            }
            
            # Making the plots
            
            all_ggplots[[k]] <- p <- ggplot(as.data.frame(post_data), aes_string(x = param_all[k]))+
              geom_histogram(aes(y = stat(count)), bins = 40, 
                             color = "black", fill = "grey") +
              #scale_y_continuous(labels = scales::percent) + 
              
              xlab(xlab)+
              ylab(y_lab) +
              
              scale_x_continuous(limits = x_lim, oob = scales::oob_keep) + 
              ylim(c(0, dim(post_data)[1])) + 
              geom_vline(xintercept = full_data[1, param_all[k]], col = "red", linetype="dotted") + #observed
              geom_vline(xintercept = quantile(post_data[, param_all[k]], 
                                               probs = 0.05, 
                                               na.rm = TRUE),
                         col="blue") + 
              geom_vline(xintercept = quantile(post_data[, param_all[k]], 
                                               probs = 0.95, 
                                               na.rm = TRUE),
                         col="blue") +
              theme_bw()
          }
          
          
          
          if(return_grid){
            return(all_ggplots)
          }else{
            grid.arrange(all_ggplots[[1]],
                         all_ggplots[[2]],
                         all_ggplots[[3]],
                         all_ggplots[[4]],
                         all_ggplots[[5]],
                         ncol = 5,
                         top = paste(z, " Adjusted Posterior", "--", 
                                     quant_x*100, "th percentile \n", dist_k, "\n (",
                                     g_union_check, " Iterations)",sep = ""))
            
            return(all_ggplots)
          }
        }
      }
      
      if(return_grid){ # On quantile and distance at a time here 
        
        return_plots <- post_plots(full_data = final,
                                   dist_k = dist_types,
                                   quant_x = quant,
                                   prior)
        
        return_plots[["contribution"]] <- contribution
        
        return(return_plots)
        
      }else{
        return(post_plots(full_data = final,
                          dist_k = dist_types,
                          quant_x = quant,
                          prior = prior))
      }
      
    }
    
    
    setwd("/Volumes/SanDisk/Research/JP/HIV/Data/ref_table")
    file_name <- paste(z, "ABC_ref_table.csv", sep = "_")
    
    write.csv(final, file = file_name)

  }
}


# This function analyses the discrepancy when recovering 
# parameters as a function of the lag component.
ABC_lag_plots <- function(request,
                          population,
                          param_all,
                          dist_types,
                          quant,
                          g_union_check,
                          stat_table = TRUE,
                          root = TRUE,
                          euclidean,
                          bins = 15){ # bins isnt actually needed 
  
  
  for(z in population){
    #keep
    if(euclidean){
      lag_dis <- matrix(NA, 0, 4) # stores discrepancy
    }else{
      # Be mindful we stopped using the standard errors so that is not included here but is above
      lag_dis <- matrix(NA, 0, 8) # stores discrepancy
    }
    for(group in c(1, 2)){
      setwd("/Volumes/SanDisk/Research/JP/HIV/Data/ref_table")
      # Generating parameter matrix
      file <- paste(z, "param_table.csv", sep = "_")
      param_matrix <- read.csv(file)
      
      
      # Reading in reference table
      
      setwd("/Volumes/SanDisk/Research/JP/HIV/Data/ref_table")
      file_name <- paste(z, "ABC_stat_table.csv", sep = "_")
      
      reg_sum_data <- read.csv(file_name)
      
      reg_sum_data <- reg_sum_data[, request, drop = FALSE] # removing extra column
      
      if(group == 1){
        reg_adj = FALSE
      }else{
        if(group == 2){
          reg_adj = TRUE
        }
      }
      
      # all lag components
      lag_c <- sort(unique(param_matrix[which(param_matrix[, "g_union"] == g_union_check),
                                        "lag"]))
      
      for(lag_i in 1:length(lag_c)){ # starting at two becuase there is no lag 0 minus obs data
        
        final <- cbind(param_matrix, reg_sum_data)
        final <- final[which(final[, "lag"] == lag_c[lag_i]), ]
        
        # Sweden
        
        observed_params_n <- sample(1:dim(final)[1], 500)
        
        observed_params <- final[observed_params_n, ]
        simulated_ref <- final[-observed_params_n, ]
        
        nn <- sample(1:dim(simulated_ref)[1], 10000)
        simulated_ref <- simulated_ref[nn, ]
        
        
        # Calculates average discrepancy over all parameters
        avg_dic <- multi_heat(observed_data = observed_params,
                              simulated_data = simulated_ref,
                              euclidean = euclidean,
                              param_all = param_all,
                              stat_table = stat_table,
                              dist_types = dist_types,
                              g_union_check = g_union_check,
                              lag_m = lag_c[lag_i],
                              z = z,
                              bins = bins,
                              plots = FALSE,
                              prior = FALSE,
                              root = root,
                              reg_adj = reg_adj)
        
        if(euclidean){
          temp <- c(lag_c[lag_i], mean(avg_dic), sd(avg_dic)/sqrt(length(avg_dic)), group)
          lag_dis <- rbind(lag_dis, temp)
        }else{
          temp <- c(lag_c[lag_i], apply(avg_dic, 2, mean), group)
          lag_dis <- rbind(lag_dis, temp)
        }

        print(lag_dis)
        #print(lag_c[1:lag_i])
        
      }
      
    }
    
    if(euclidean){
      colnames(lag_dis) <- c("lag", "mean", "s.e", "group")
    }else{
      colnames(lag_dis) <- c("lag", paste("mean", param_all[1:5], sep = "_"), "all", "group")
    }
    
    setwd("/Volumes/SanDisk/Research/JP/HIV/Data/figures")
    if(root){
      write.csv(lag_dis, paste("Root_Euclidean_", euclidean, "_RMSE_data.csv", sep = ""), row.names = FALSE)
      
    }else{
      write.csv(lag_dis, paste("Euclidean_", euclidean, "_RMSE_data.csv", sep = ""), row.names = FALSE)
    }

    final <- cbind(param_matrix, reg_sum_data)
    final <- final[which(final[, "lag"] == -1), ] # lag should not matter because we are using full prior
    
    # # Sweden
    #observed_params_n <- sample(which(final[, "sigma"] < 50 & final[, "sigma"] >= 40), 
    #                            size = 100, 
    #                            replace = FALSE)
    # 
    observed_params_n <- sample(1:dim(final)[1], 500)
    # 
    observed_params <- final[observed_params_n, ]
    simulated_ref <- final[-observed_params_n, ]
    # # Generating discrepancy from the prior
    prior_point <- multi_heat(observed_data = observed_params,
                              simulated_data = simulated_ref,
                              euclidean = TRUE,
                              param_all = c(param_all, "mu", "lag", "g_union"),
                              stat_table = stat_table,
                              dist_types = dist_types,
                              g_union_check = g_union_check,
                              lag_m = -1, # one time point dependent 
                              z = z,
                              bins = bins,
                              plots = FALSE,
                              prior = TRUE,
                              root = TRUE,
                              reg_adj = TRUE) # but we do prior as true
    
    if(euclidean){
      colnames(lag_dis) <- c("lag", "mean", "s.e", "group")
    }else{
      colnames(lag_dis) <- c("lag", paste("mean", param_all[1:5], sep = "_"), "all", "group")
    }
    lag_data <- as.data.frame(lag_dis)
    lag_data[, "group"] <- as.factor(lag_data[, "group"] + 1) 
    
    
    if(euclidean){
      p <- ggplot(lag_data, aes(x=lag, y=mean, group = group, colour = group)) + 
        geom_point(size = 5) +
        ylab("Average Discrepancy") +
        xlab("Iteration Lag") +
        ggtitle(paste("Stockholm", " Lag Accuracy Plots (",
                      g_union_check, " iteration(s)) \n",
                      paste(request, collapse =", ") ,sep = "")) + 
        geom_smooth(method = "loess", size = 1.5) + 
        theme(axis.title.x = element_text(size = 20),
              axis.title.y = element_text(size = 20)) + 
        theme_bw() +
        geom_hline(yintercept= mean(prior_point), linetype="dashed", color = "red") + 
        scale_color_manual(labels = c("noAdj", "Adj"),
                           values = c("blue", "red"))
    }else{
      mdf <- melt(lag_data,id.vars= c("lag", "group"))
      
      p <- ggplot(mdf, aes( x=lag, y=value, colour=variable, group= variable)) + 
        geom_point(size = 5)  + 
        geom_smooth(method = "loess", size = 0.9) +
        ylab("Average Discrepancy") +
        xlab("Observation Lag") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, size = 18),
              axis.text.y = element_text(size = 18),
              text = element_text(size = 18)) +
        theme(legend.position="top") +
        #geom_hline(yintercept= prior_point, linetype="dashed", color = "red") + 
        scale_color_manual(labels = expression(rho,sigma, omega[0], omega[1],
                                               "Total"), name = NULL,
                           values = c("blue", "red", "green", "grey","black")) 
    }
    
    return(list(p, lag_data, request))
    
  }
  
}


ABC_lag_plots(request,
              population,
              param_all,
              dist_types,
              quant,
              g_union_check,
              stat_table = FALSE,
              root = TRUE,
              euclidean = FALSE,
              bins)

ABC_lag_plots(request,
              population,
              param_all,
              dist_types,
              quant,
              g_union_check,
              stat_table = FALSE,
              root = TRUE,
              euclidean = TRUE,
              bins)

# We use this function to make plots if our base files are 
# from csv
plots_csv <- function(iter_1, 
                      request,
                      population,
                      param_all,
                      d_c){
  # iter_1:
    # The number of iterations to consider
  # request:
    # Summary statistics we are requesting
    # num_edges: number of edges
    # edges_m: proportion of edges maintained between time points
    # num_tri: number of triangles in the graph
    # quant_25: 25th degree distribution percentile
    # quant_50: 50th degree distribution percentile
    # quant_75: 75th degree distribution percentile
    # quant_max: maximum of the degree distribution 
    # quant_mean: mean of the degree distribution 
    # median_local_cc: median of the local clustering coefficient
    # mean_local_cc: mean of the local clustering coefficient 
    # num_con_comp: number of components with 2 or more nodes
    # size_LCC: size of the largest connected component
    # num_singletons: number of nodes with degree 0
    # num_exclusive_pairs: number of monogamous pairs
  # population:
    # Country/City we are looking at
      # Sweden
      # Amsterdam
  # param_all
    # parameters to be considered (each model will be different)
  # d_c: "discrete_free" or "discrete_fixed"
  all_graphs <- list()
  for(z in population){
    for(request_i in request){
      for(i in param_all){
        name_i <- i
        print(c(z, request_i, i))
        wd_all <- file.path("/Volumes/SanDisk/Research/JP/HIV/Results/Sum_Discovery",
                            z,
                            "Time_Varying",
                            i,
                           d_c)
        setwd(wd_all)
        
        # Reads in all files in the directory
        all_files <- list.files()
        
        # Stores all of the iteration varying results independent of time
        # varying
        iter_removed <- character_remove(all_files)
        dim_all <- matrix(NA, iter_1, length(iter_removed)) # change column
        colnames(dim_all) <- processor(iter_removed)[[1]]
        iter_removed <- iter_removed[processor(iter_removed)[[2]]] # ordering the names by value size
        
        for(k in 1:dim(dim_all)[2]){ # iterates through base parameter
          print(k/dim(dim_all)[2])
          # and parameter value combos
          
          # Additional "_" is added to seperate 1 from 11, 1.1 etc
          all_files_k <- grep(paste(iter_removed[k], "_", sep = ""), 
                              all_files, 
                              value = TRUE)
          
          dim_iter <- rep(NA, iter_1)
          for(x in 1:I(length(all_files_k)/2)){# iterates through the iterations.
            # Divide by 2 to account for time varying files
            
            all_files_k_iter <- grep(paste("^",x, "_.*$", sep = ""), 
                                     all_files_k,
                                     value = TRUE)
            
            if(length(all_files_k_iter) == 0){
              next
            }
            
            G <- list() # stores time dependent graphs
            edge_type_data <- list() # stores edge type data
            for(y in 0:1){ # iterates through time points
              #print(c(k, x, y))
              vertex_search <- gsub("edgeList", 
                                    "vertexList",
                                    all_files_k_iter[y + 1])
              
              vertex_string_procc <- paste("^", vertex_search, "$", sep = "")
              
              vertex_file <- grep(vertex_string_procc,
                                  all_files,
                                  value = TRUE)
              
              vertex_data <- read.csv(vertex_file)
              edge_data <- read.csv(all_files_k_iter[y + 1])
              
              G[[y + 1]] <- make_empty_graph(n = dim(vertex_data)[1],
                                             directed = FALSE)
              
              V(G[[y + 1]])$name <- vertex_data[, 1]
              
              G[[y + 1]] <- add_edges(graph = G[[y + 1]],
                                      edges = c(t(edge_data)))
              
              if(y == 0){
                edge_type_file <- gsub("vertexList",
                                       "steady_casual_tracker",
                                       vertex_file)
                edge_type_data[[y+1]] <- read.csv(edge_type_file)
              }
              
            }
            
            # Calculating the summary statistic of interest
            dim_iter[x] <- summ_stat_calc(G = G,
                                          request = request_i,
                                          steady_casual_files = edge_type_data,
                                          lag = 15)
            
          }
          
          dim_all[, k] <- dim_iter
        }
        
        
        dim_all <- as.data.frame(dim_all)
        setwd("/Volumes/SanDisk/Research/JP/HIV/Data/figures")
        write.csv(dim_all, paste(i, request_i , d_c, "mapping_data.csv", sep = "_"), row.names = FALSE)
        
        # Creating ggplot
        #title_1 <- paste(z, "Summary Statistic Plot")
        title_1 <- NULL
        if(request_i == "num_edges"){
          ylab_1 <- "# of Edges"
        }else{
          if(request_i == "edges_m"){
            ylab_1 <- "Proportion of Edges Maintained in Final Time Point"
          }else{
            if(request_i == "num_tri"){
              ylab_1 <- "# of Triangles"
            }else{
              if(request_i == "quant_25"){
                ylab_1 <- "25th Degree Distribution Percentile"
              }else{
                if(request_i == "quant_50"){
                  ylab_1 <- "50th Degree Distribution Percentile"
                }else{
                  if(request_i == "quant_75"){
                    ylab_1 <- "75th Degree Distribution Percentile"
                  }else{
                    if(request_i == "quant_max"){
                      ylab_1 <- "Max of Degree Distribution"
                    }else{
                      if(request_i == "quant_mean"){
                        ylab_1 <- "Mean of Degree Distribution"
                      }else{
                        if(request_i == "median_local_cc"){
                          ylab_1 <- "Median of the Local Clustering Coefficient"
                        }else{
                          if(request_i == "mean_local_cc"){
                            ylab_1 <- "Mean of the Local Clustering Coefficient"
                          }else{
                            if(request_i == "num_con_comp"){
                              ylab_1 <- "Number of Connected Components (2 or More Nodes)"
                            }else{
                              if(request_i == "size_LCC"){
                                ylab_1 <- "Size of Largest Connected Component"
                              }else{
                                if(request_i == "num_singletons"){
                                  ylab_1 <- "num_singletons"
                                }else{
                                  if(request_i == "num_exclusive_pairs"){
                                    ylab_1 <- "Number of Monogamous Relationships"
                                  }else{
                                    if(request_i == "concurrency"){
                                      ylab_1 <- "Proportion of Individuals With Multiple Partners"
                                    }else{
                                      if(request_i == "verticies_m"){
                                        ylab_1 <- "Proportion of Verticies Maintained in Final Time Point"
                                      }else{
                                        if(request_i == "prolonged_single"){
                                          ylab_1 <- "Proportion of individuals who stay single"
                                        }else{
                                          if(request_i == "entering"){
                                            ylab_1 <- "Proportion of New Individuals in the Population"
                                          }else{
                                            ylab_1 <- request_i
                                          }
                                        }
                                      }
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
        
        
        
        if(i == "p_cas"){
          i <-  "Probability a Casual Relationship is Maintained"
        }else{
          if(i == "k_max"){
            i <- "Maximum Allowable Degree"
          }else{
            if(i == "w_0"){
              i <- "w_0"
            }else{
              if(i == "w_1"){
                i <- "w_1"
              }else{
                if(i == "migration"){
                  i <- "Rate Individuals Age Into the Population"
                }else{
                  if(i == "p_0"){
                    i <- "Probability a Node has Degree 0"
                  }else{
                    i <- parse(text = i)
                  }
                }
              }
            }
          }
        }
        
        
        p1 <- ggplot(stack(dim_all), aes(x = ind, y = values))+
          theme_bw() +
          geom_boxplot() + 
          theme(axis.text.x = element_text(angle = 90, size = 6),
                axis.text.y = element_text(size = 6)) +
          xlab(i) +
          ylab(ylab_1) + 
          labs(title = title_1) 
        
        plot(p1)
        
        all_graphs[[name_i]][[request_i]] <- p1
      }
    }
  }
  
  grid.arrange(grobs = c(all_graphs[[1]], # not doing mu
                         all_graphs[[2]], 
                         all_graphs[[3]],
                         all_graphs[[4]]), nrow = 4, ncol = length(request))
  
  return(all_graphs)
}


# Mapping functions
plots_csv(iter_1 = 100, 
          request, 
          population, 
          param_all[2:5],
          d_c = "discrete_fixed")

plots_csv(iter_1 = 100, 
          request, 
          population, 
          param_all[2:5],
          d_c = "discrete_free")
  

