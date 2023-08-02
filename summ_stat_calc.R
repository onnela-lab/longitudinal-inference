################################################################################################
# This function calculates different summary statistics depending on the input.
# t_varying should be TRUE for any graph inputed in this algorithm.

# Date started: 10/14/2021
# Last edited: 10/14/2021
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
  # num_edges: number of edges
  # edges_m: proportion of edges maintained between time points
  # num_tri: number of triangles in the graph
  # quant_25: 25th degree distribution quantile
  # quant_50: 50th degree distribution quantile
  # quant_75: 75th degree distribution quantile
  # quant_max: maximum of the degree distribution 
  # quant_mean: mean of the degree distribution 
  # median_local_cc: median of the local clustering coefficient
  # mean_local_cc: mean of the local clustering coefficient 
  # num_con_comp: number of components with 2 or more nodes
  # size_LCC: size of the largest connected component
  # num_singletons: number of nodes with degree 0
  # num_exclusive_pairs: number of monogamous pairs
  # verticies_m: proportion of first iterations in the last iteration
  # concurrency: proportion of nodes with more than one edge
  # prolonged_single: Proportion of single individuals still single
  # entering: propotion of new individuals in the population 
  # mean_distance: Average length of all of the shortest paths
  # single_partner_ratio: number of single individuals divided number of edges
  # part_length_single_ratio: Proportion of relationships maintained/number of single individuals.
  # seperate_migrate_ratio: seperation percentage divided by number of individulas who left
  # inverse_sep_minus_migrate: inverse seperation minus twice migration number
  # steady_m: proportion of maintained steady relationships maintained
  # casual_m: proportion of maintained casual relationships maintained
  # two_concurrency: mean percentage of people in steady relationships also in casual ones
  # new_edges: what percentage of present steady edges are new
  # clique_mean: average size of all of the cliques
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
  
  # one day casual iterations makes sure we always get the last iteration
  last_iteration <- steady_casual_files[[1]][2, "iteration"] # couldve been one *shrugs*
  #last_iteration <- max(steady_casual_files[[1]][, "start"])
  
  # lag -1 implies same time point. 
  lag <- lag + 1
  
  if(request == "num_edges"){
    sum_stat_0 <- dim(E_0)[1]
    sum_stat_1 <- dim(E_1)[1]
    #print(E_0)
    #print(dim(E_0)[1])
    
    #print(E_1)
    #print(dim(E_1)[1])
    
    sum_stat <- (sum_stat_0 + sum_stat_1)/2
    
    #print(sum_stat)
  }else{
    if(request == "edges_m"){
      c_edge <- common_edges(E_0, E_1)
      sum_stat <- (dim(c_edge)[1]/dim(E_0)[1])*100
    }else{
      if(request == "num_tri"){
        sum_stat_0 <- length(triangles(g_0))/3
        sum_stat_1 <- length(triangles(g_1))/3
        
        sum_stat <- (sum_stat_0 + sum_stat_1)/2
      }else{
        if(request == "quant_25"){
          sum_stat_0 <- as.numeric(summary(degree(g_0))["1st Qu."])
          sum_stat_1 <- as.numeric(summary(degree(g_1))["1st Qu."])
          
          sum_stat <- (sum_stat_0 + sum_stat_1)/2
        }else{
          if(request == "quant_50"){
            sum_stat_0 <- as.numeric(summary(degree(g_0))["Median"])
            sum_stat_1 <- as.numeric(summary(degree(g_1))["Median"])
            
            sum_stat <- (sum_stat_0 + sum_stat_1)/2
          }else{
            if(request == "quant_75"){
              sum_stat_0 <- as.numeric(summary(degree(g_0))["3rd Qu."])
              sum_stat_1 <- as.numeric(summary(degree(g_1))["3rd Qu."])
              
              sum_stat <- (sum_stat_0 + sum_stat_1)/2
            }else{
              if(request == "quant_max"){
                sum_stat_0 <- as.numeric(summary(degree(g_0))["Max."])
                sum_stat_1 <- as.numeric(summary(degree(g_1))["Max."])
                
                sum_stat <- (sum_stat_0 + sum_stat_1)/2
              }else{
                if(request == "quant_mean"){
                  sum_stat_0 <- as.numeric(summary(degree(g_0))["Mean"])
                  sum_stat_1 <- as.numeric(summary(degree(g_1))["Mean"])
                  
                  sum_stat <- (sum_stat_0 + sum_stat_1)/2
                }else{
                  if(request == "median_local_cc"){
                    
                    trans_0 <- transitivity(graph = g_0,
                                            type = c("local"),
                                            vids = NULL,
                                            weights = NULL,
                                            isolates = c("NaN")
                    )
                    
                    sum_stat_0 <- median(trans_0, na.rm = TRUE)
                    
                    trans_1 <- transitivity(graph = g_1,
                                            type = c("local"),
                                            vids = NULL,
                                            weights = NULL,
                                            isolates = c("NaN")
                    )
                    
                    sum_stat_1 <- median(trans_1, na.rm = TRUE)
                    
                    
                    sum_stat <- (sum_stat_0 + sum_stat_1)/2
                  }else{
                    if(request == "mean_local_cc"){
                      
                      trans_0 <- transitivity(graph = g_0,
                                              type = c("local"),
                                              vids = NULL,
                                              weights = NULL,
                                              isolates = c("NaN")
                      )
                      
                      sum_stat_0 <- mean(trans_0, na.rm = TRUE)
                      
                      trans_1 <- transitivity(graph = g_1,
                                              type = c("local"),
                                              vids = NULL,
                                              weights = NULL,
                                              isolates = c("NaN")
                      )
                      
                      sum_stat_1 <- mean(trans_1, na.rm = TRUE)
                      
                      
                      sum_stat <- (sum_stat_0 + sum_stat_1)/2
                    }else{
                      if(request == "num_con_comp"){
                        sum_stat_0 <- length(decompose(graph = g_0, 
                                                       mode = "strong", 
                                                       min.vertices = 2))
                        
                        sum_stat_1 <- length(decompose(graph = g_1, 
                                                       mode = "strong", 
                                                       min.vertices = 2))
                        
                        sum_stat <- (sum_stat_0 + sum_stat_1)/2
                      }else{
                        if(request == "size_LCC"){
                          sum_stat_0 <- max(clusters(g_0, mode="strong")$csize)
                          sum_stat_1 <- max(clusters(g_1, mode="strong")$csize)
                          
                          sum_stat <- (sum_stat_0 + sum_stat_1)/2
                        }else{
                          if(request == "num_singletons"){
                            sum_stat_0 <- sum(degree(g_0) == 0)
                            sum_stat_1 <- sum(degree(g_1) == 0)
                            
                            sum_stat_0 <- sum_stat_0/length(V(g_0))
                            sum_stat_1 <- sum_stat_1/length(V(g_1))
                            
                            sum_stat <- (sum_stat_0 + sum_stat_1)/2
                          }else{
                            if(request == "num_exclusive_pairs"){
                              sum_stat_0 <- sum(clusters(g_0, mode="strong")$csize == 2)
                              sum_stat_1 <- sum(clusters(g_1, mode="strong")$csize == 2)
                              
                              sum_stat <- (sum_stat_0 + sum_stat_1)/2
                            }else{
                              if(request == "verticies_m"){
                                sum_stat <- mean(V(g_0)$name %in% V(g_1)$name)*100
                              }else{
                                if(request == "h_concurrency"){# greater than mean and more than one partner
                                  sum_stat_0 <- mean(degree(g_0) > mean(degree(g_0)) &
                                                       degree(g_0) > 1)*100
                                  
                                  sum_stat_1 <- mean(degree(g_1) > mean(degree(g_1)) &
                                                       degree(g_1) > 1)*100
                                  
                                  sum_stat <- (sum_stat_0 + sum_stat_1)/2
                                }else{
                                  if(request == "l_concurrency"){# greater than mean and more than one partner
                                    sum_stat_0 <- mean(degree(g_0) > 1)*100
                                    sum_stat_1 <- mean(degree(g_1) > 1)*100
                                    
                                    sum_stat <- (sum_stat_0 + sum_stat_1)/2
                                  }else{
                                    if(request == "prolonged_single"){
                                      initial_single <- names(which(degree(g_0) == 0))
                                      end_single <- names(which(degree(g_1) == 0))
                                      
                                      if(length(initial_single) == 0|length(end_single) == 0){
                                        still_single <- 0
                                      }else{
                                        still_single <- mean(initial_single %in% end_single)*100
                                      }
                                      
                                      sum_stat <- still_single
                                    }else{
                                      if(request == "entering"){
                                        sum_stat <- mean(!V(g_1)$name %in% V(g_0)$name)*100
                                      }else{
                                        if(request == "mean_distance"){
                                          sum_stat_0 <- mean_distance(g_0, 
                                                                      directed = TRUE)
                                          
                                          sum_stat_1 <- mean_distance(g_1, 
                                                                      directed = TRUE)
                                          
                                          sum_stat <- (sum_stat_0 + sum_stat_1)/2
                                        }else{
                                          if(request == "single_partner_ratio"){
                                            sum_stat_0 <- sum(degree(g_0) == 0)/dim(E_0)[1]
                                            sum_stat_1 <- sum(degree(g_1) == 0)/dim(E_1)[1]
                                            
                                            sum_stat <- (sum_stat_0 + sum_stat_1)/2
                                          }else{
                                            if(request == "part_length_single_ratio"){
                                              c_edge <- common_edges(E_0, E_1)
                                              part_maintained <- (100 -(dim(c_edge)[1]/dim(E_0)[1])*100)
                                              
                                              if(sum(degree(g_1) == 0) > 0){
                                                sum_stat <- part_maintained/sum(degree(g_1) == 0)
                                              }else{
                                                sum_stat <- NA
                                              }
                                            }else{
                                              if(request == "seperate_migrate_ratio"){
                                                left <- sum(!V(g_0)$name %in% V(g_1)$name)*100
                                                
                                                c_edge <- common_edges(E_0, E_1)
                                                
                                                if(dim(E_0)[1] == 0){
                                                  seperated <- 0
                                                }else{
                                                  seperated <- 100 - (dim(c_edge)[1]/dim(E_0)[1])*100
                                                }
                                                
                                                if(seperated > 0){
                                                  sum_stat <- left/seperated
                                                }else{
                                                  sum_stat <- NA
                                                }
                                              }else{
                                                if(request == "inverse_sep_minus_migrate"){
                                                  c_edge <- common_edges(E_0, E_1)
                                                  if(dim(E_0)[1] == 0){
                                                    seperated <- 0
                                                  }else{
                                                    seperated <- 100 - (dim(c_edge)[1]/dim(E_0)[1])*100
                                                    left <- sum(!V(g_0)$name %in% V(g_1)$name)*100
                                                  }
                                                  
                                                  
                                                  if(seperated > 0){
                                                    sum_stat <- 1/seperated - 2*left
                                                  }else{
                                                    sum_stat <- NA
                                                  }
                                                  
                                                }else{
                                                  if(request %in% c("negative_change", "positive_change", "no_change")){
                                                    name_0 <- V(g_0)$name[which(V(g_0)$name %in% V(g_1)$name)]
                                                    name_all <- V(g_1)$name[V(g_1)$name %in% name_0]
                                                    
                                                    changed_degree <- degree(g_1)[name_all] - degree(g_0)[name_all]
                                                    
                                                    if(request == "negative_change"){
                                                      sum_stat <- sum(changed_degree < 0)/length(name_all)
                                                    }else{
                                                      if(request == "positive_change"){
                                                        sum_stat <- sum(changed_degree > 0)/length(name_all)
                                                      }else{
                                                        if(request == "no_change"){
                                                          sum_stat <- sum(changed_degree == 0)/length(name_all)
                                                        }else{
                                                          if(request == "single_num_by_single_length"){
                                                            initial_single <- names(which(degree(g_0) == 0))
                                                            end_single <- names(which(degree(g_1) == 0))
                                                            
                                                            if(length(initial_single) == 0|length(end_single) == 0){
                                                              still_single <- 0
                                                            }else{
                                                              still_single <- mean(initial_single %in% end_single)
                                                            }
                                                            
                                                            single_num <- mean(degree(g_1) == 0)
                                                            
                                                            if(single_num != 0 & still_single != 0){
                                                              sum_stat <- 1/(single_num*still_single)
                                                            }else{
                                                              sum_stat <- 0
                                                            }
                                                          }
                                                        }
                                                      }
                                                    }
                                                  }else{
                                                    if(request == "steady_end_len"){
                                                      if(is.null(steady_casual_files)){
                                                        sum_stat <- NA
                                                        
                                                      }else{
                                                        initial_steady <- steady_casual_files[[1]][which(steady_casual_files[[1]][, "type"]
                                                                                                         == "steady" | 
                                                                                                           steady_casual_files[[1]][, "type"]
                                                                                                         == "s"), ,drop = FALSE]
                                                        
                                                        # which relationships were ongoing during duration of study and ends
                                                        # during duration of the study
                                                        nn <- which(initial_steady[, "end"] >= last_iteration - lag - 3 & # added 3 iterations to first window
                                                                      initial_steady[, "end"] <= last_iteration)
                                                        #print(length(nn))
                                                        
                                                        if(length(nn) > 0){
                                                          duration <- initial_steady[nn, "end"] - initial_steady[nn, "start"]
                                                          sum_stat <- mean(duration, na.rm = TRUE)
                                                        }else{
                                                          sum_stat <- NA
                                                        }  
                                                      }
                                                    }else{
                                                      if(request == "casual_len"){
                                                        if(is.null(steady_casual_files)){
                                                          sum_stat <- NA
                                                          
                                                        }else{
                                                          
                                                          initial_steady <- steady_casual_files[[1]][which(steady_casual_files[[1]][, "type"]
                                                                                                           == "casual" | 
                                                                                                             steady_casual_files[[1]][, "type"]
                                                                                                           == "c"), ,drop = FALSE]
                                                          
                                                          # which relationships were ongoing during duration of study and ends
                                                          # during duration of the study
                                                          nn <- which(initial_steady[, "end"] >= last_iteration - lag &
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
                                                          if(request %in% c("random_1", "random_2", "random_3",
                                                                            "random_4", "random_5", "random_6", 
                                                                            "random_7", "random_8",
                                                                            "random_9", "random_10")){
                                                            sum_stat <- runif(1, 0, 5)
                                                          }else{
                                                            if(request %in% c("random_normal")){
                                                              sum_stat <- rnorm(1)
                                                            }else{
                                                              if(request == "new_edges"){
                                                                initial_steady <- steady_casual_files[[1]][which(steady_casual_files[[1]][, "type"]
                                                                                                                 == "steady" | 
                                                                                                                   steady_casual_files[[1]][, "type"]
                                                                                                                 == "s"), ,drop = FALSE]
                                                                
                                                                # What proportion of final graph is new edges
                                                                nn_0 <- which(initial_steady[, "start"] <= last_iteration - lag & 
                                                                                initial_steady[, "end"] > last_iteration - lag)
                                                                
                                                                nn_1 <- which(initial_steady[, "end"] > last_iteration)
                                                                
                                                                c_edge <- common_edges(initial_steady[nn_0, c("e1", "e2")], 
                                                                                       initial_steady[nn_1, c("e1", "e2")])
                                                                
                                                                if(length(nn_1) == 0){
                                                                  sum_stat <- 0
                                                                }else{ # common edges represent edges that are the same between graphs
                                                                  sum_stat <- 100 - dim(c_edge)[1]/length(nn_1)*100
                                                                }
                                                                
                                                              }else{
                                                                if(request == "steady_num"){
                                                                  end_steady_0 <- steady_casual_files[[1]][which(steady_casual_files[[1]][, "type"]
                                                                                                                 == "steady" | 
                                                                                                                   steady_casual_files[[1]][, "type"]
                                                                                                                 == "s"), ,drop = FALSE]
                                                                  
                                                                  #print(end_steady_0)
                                                                  #print(end_steady_0[, "start"])
                                                                  #print(end_steady_0[, "end"])
                                                                  #print(c(last_iteration, lag))
                                                                  
                                                                  nn_0 <- which(end_steady_0[, "start"] <= last_iteration - lag & 
                                                                                  end_steady_0[, "end"] > last_iteration - lag)
                                                                  #print(end_steady_0[nn_0, ])
                                                                  #print(length(nn_0))
                                                                  
                                                                  sum_stat_0 <- length(nn_0)
                                                                  
                                                                  #print(sum_stat_0)
                                                                  
                                                                  nn_1 <- which(end_steady_0[, "end"] > last_iteration)
                                                                  #print(end_steady_0[nn_1, ])
                                                                  
                                                                  sum_stat_1 <- length(nn_1)
                                                                  
                                                                  sum_stat_0 <- sum_stat_0/dim(E_0)[1]
                                                                  sum_stat_1 <- sum_stat_1/dim(E_1)[1]
                                                                  
                                                                  
                                                                  sum_stat <- (sum_stat_0 + sum_stat_1)/2
                                                                }else{
                                                                  if(request == "casual_num"){
                                                                    end_casual_0 <- steady_casual_files[[1]][which(steady_casual_files[[1]][, "type"]
                                                                                                                   == "casual" | 
                                                                                                                     steady_casual_files[[1]][, "type"]
                                                                                                                   == "c"), ,drop = FALSE]
                                                                    
                                                                    #print(end_casual_0)
                                                                    #print(end_casual_0[, "end"])
                                                                    
                                                                    
                                                                    nn_0 <- which(end_casual_0[, "start"] == last_iteration - lag)
                                                                    #print(length(nn_0))
                                                                    
                                                                    sum_stat_0 <- length(nn_0)
                                                                    #print(end_casual_0[nn_0, ])
                                                                    
                                                                    
                                                                    nn_1 <- which(end_casual_0[, "start"] == last_iteration)
                                                                    
                                                                    sum_stat_1 <- length(nn_1)
                                                                    
                                                                    sum_stat_0 <- sum_stat_0/dim(E_0)[1]
                                                                    sum_stat_1 <- sum_stat_1/dim(E_1)[1]
                                                                    
                                                                    
                                                                    
                                                                    sum_stat <- (sum_stat_0 + sum_stat_1)/2
                                                                  }else{
                                                                    if(request == "steady_casual_ratio"){
                                                                      end_steady_0 <- steady_casual_files[[1]][which(steady_casual_files[[1]][, 3]
                                                                                                                     == "steady"|
                                                                                                                       steady_casual_files[[1]][, 3]
                                                                                                                     == "s"), 
                                                                                                               1:2, ,drop = FALSE]
                                                                      
                                                                      st_0 <- c(end_steady_0[, 1], end_steady_0[, 2])
                                                                      
                                                                      end_casual_0 <- steady_casual_files[[1]][which(steady_casual_files[[1]][, 3]
                                                                                                                     == "casual"|
                                                                                                                       steady_casual_files[[1]][, 3]
                                                                                                                     == "c"),
                                                                                                               1:2, ,drop = FALSE]
                                                                      
                                                                      if(dim(end_casual_0)[1] != 0){
                                                                        sum_stat_0 <- dim(end_steady_0)[1]/dim(end_casual_0)[1]
                                                                      }else{
                                                                        sum_stat_0 <-  dim(end_steady_0)[1]
                                                                      }
                                                                      
                                                                      
                                                                      end_steady_1 <- steady_casual_files[[2]][which(steady_casual_files[[2]][, 3]
                                                                                                                     == "steady"|
                                                                                                                       steady_casual_files[[2]][, 3]
                                                                                                                     == "s"), 
                                                                                                               1:2, ,drop = FALSE]
                                                                      
                                                                      st_1 <- c(end_steady_1[, 1], end_steady_1[, 2])
                                                                      
                                                                      end_casual_1 <- steady_casual_files[[2]][which(steady_casual_files[[2]][, 3]
                                                                                                                     == "casual"|
                                                                                                                       steady_casual_files[[2]][, 3]
                                                                                                                     == "c"),
                                                                                                               1:2, ,drop = FALSE]
                                                                      
                                                                      if(dim(end_casual_1)[1] != 0){
                                                                        sum_stat_1 <- dim(end_steady_1)[1]/dim(end_casual_1)[1]
                                                                      }else{
                                                                        sum_stat_1 <-  dim(end_steady_1)[1]
                                                                      }
                                                                      
                                                                      sum_stat <- (sum_stat_0 + sum_stat_1)/2
                                                                    }else{
                                                                      if(request == "clique_mean"){
                                                                        
                                                                        clique_0 <- cliques(g_0) 
                                                                        clique_1 <- cliques(g_1) 
                                                                        
                                                                        clique_size_0 <- NULL
                                                                        clique_size_1 <- NULL
                                                                        
                                                                        for(i in 1:length(clique_0)){
                                                                          clique_size_0[i] <- length(clique_0[[i]])
                                                                        }
                                                                        
                                                                        for(i in 1:length(clique_1)){
                                                                          clique_size_1[i] <- length(clique_1[[i]])
                                                                        }
                                                                        
                                                                        sum_stat_0 <- mean(clique_size_0)
                                                                        sum_stat_1 <- mean(clique_size_1)
                                                                        
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
                                                                            #print(length(nn))
                                                                            
                                                                            if(length(nn) > 0){
                                                                              duration <- initial_steady[nn, "end"] - initial_steady[nn, "start"]
                                                                              sum_stat <- mean(duration, na.rm = TRUE)
                                                                            }else{
                                                                              sum_stat <- NA
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
    }
  }
  
  return(sum_stat) 
}


#summ_stat_calc(G, "edges_m")
