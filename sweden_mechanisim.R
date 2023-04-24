########################################################################################################################
# This script codes the mechanistic model for the Sweden population. HIV transmission is not considered in this 
# script.

# All rates are considered to be the parameter of geometric distribution, except the constant rate
# at which people leave the population.

# We only consider the discrete distribution, and we do not add nor remove nodes due to
# fixing mu = 0 for a closed population.
########################################################################################################################
source("/home/ot25/Research/JP/HIV/Scripts/edge_index.R")

Sweden_MSM <- function(num_nodes,
                       mu,
                       rho,
                       sigma,
                       w_0,
                       w_1,
                       iter,
                       g_union,
                       distribution,
                       lag = 0){
  
  library(igraph)
  
  ## This section contains all functions used in the loop
  
  # General edge remove function
  Edge_Remove <- function(G, remove){
    #G: Graph to be modified
    #remove: edge sequence to be removed
    n_e <- edge_index(G = G,
                      edge_sequence = c(t(remove)))
    
    G <- delete_edges(graph = G,
                      edges = n_e)
    
    return(G)
  }
  
  # Casual/Steady edge remove function
  Casual_Steady_Remove <- function(G, dynamic_edge_tracker, i, sigma, distribution, type){
    # G: Graph to be modified
    # dynamic_edge_tracker: Object tracking all edge information
    # i: Present iteration
    # sigma: Edge removal probability
    # distribution: Continuous or discrete
    # type: Type of edges we will be removing
    
    if(type == "c"){ # Removal process for casual edges
      remove <- dynamic_edge_tracker[which(as.numeric(dynamic_edge_tracker[, "end"]) == i  &
                                             dynamic_edge_tracker[, "type"] == "c"), c("e1", "e2"), 
                                     drop = FALSE]
    }else{
      if(type == "s"){ # Removal process for steady edges
        if(distribution == "continuous"){
          n_sep_2 <- which(as.numeric(dynamic_edge_tracker[, "end"]) == i & 
                             dynamic_edge_tracker[, "type"] == "s")
          remove <- dynamic_edge_tracker[n_sep_2, c("e1", "e2"), drop = FALSE]
          
        }else{
          if(distribution == "discrete"){ # We only consider the discrete distribution
            n_current_edges <- which(dynamic_edge_tracker[, "end"] == "Inf" & # edges still present
                                       dynamic_edge_tracker[, "type"] == "s")
            
            if(length(n_current_edges) > 0){
              current_edges <- dynamic_edge_tracker[n_current_edges, c("e1", "e2"), drop = FALSE]
              n_r_e <- c(1:dim(current_edges)[1])[which(runif(dim(current_edges)[1]) < sigma)]
              
              remove <- current_edges[n_r_e, , drop = FALSE]
            }else{ # We dont remove edges unless we get at least one coin flip saying to
              remove <- matrix(NA, 0, 1)
            }
          }
        }
      }
    }
    
    if(dim(remove)[1] > 0){
      G <- Edge_Remove(G, remove)
      if(distribution == "discrete" & type == "s"){ # Tracking when steady edges end
        dynamic_edge_tracker[n_current_edges[n_r_e], "end"] <- i
      }
    }
    
    
    return(list(G, dynamic_edge_tracker))
  }
  
  # Node remove
  Node_Remove <- function(G, dynamic_edge_tracker, exit_time, i, mu, distribution){
    # G: Graph to be modified
    # dynamic_edge_tracker: Object tracking all edge information
    # exit_time: Each nodes exit time (not useful in discrete case)
    # i: Present iteration
    # mu: node removal probability
    # distribution: Continuous or discrete
    
    if(distribution == "continuous"){
      # Deleting associated edges with those nodes being removed
      n_r <- names(which(exit_time == i))
      
    }else{
      if(distribution == "discrete"){
        # Deleted edges are discretely decided
        n_r <- V(G)$name[which(runif(length(V(G))) < mu)]
      }
    }
    
    if(length(n_r) > 0){
      for(k in 1:length(n_r)){ # this is done differently because we are identifying based on nodes not edges
        n_sep <- edge_index(G = G,
                            edge_sequence = n_r[k])
        
        if(length(n_sep) > 0){ # Deleting associated edges.
          G <- delete_edges(graph = G,
                            edges = n_sep)
          
          n_sep_r <- which((dynamic_edge_tracker[, "e1"] %in% n_r[k]|
                              dynamic_edge_tracker[, "e2"] %in% n_r[k]) &
                             as.numeric(dynamic_edge_tracker[, "end"]) > i - 1 & # edge still present
                             dynamic_edge_tracker[, "type"] == "s") 
          
          # The new end time is i
          dynamic_edge_tracker[n_sep_r, "end"] <- i # steady edges are ended by migration
          
        }
        
      }
      
      # deleting nodes
      G <- delete_vertices(G, v = n_r)
    }
    
    return(list(G, dynamic_edge_tracker))
    
  }
  
  # Add Nodes
  Node_Add <- function(G, dynamic_edge_tracker, exit_time, i, mu, distribution){
    # G: Graph to be modified
    # dynamic_edge_tracker: Object tracking all edge information
    # exit_time: Each nodes exit time (not useful in discrete case)
    # i: Present iteration
    # mu: node removal probability
    # distribution: Continuous or discrete
    
    # Deciding how many nodes to add. 4.5 is either 4 or 5.
    dec <- (num_nodes * mu) - floor(num_nodes * mu)
    
    if(runif(1) < dec){
      nodes_add <- ceiling(num_nodes * mu)
    }else{
      nodes_add <- floor(num_nodes * mu)
    }
    
    if(nodes_add == 0){
      return(list(G, dynamic_edge_tracker, exit_time))
    }
    
    # Adding nodes to the graph
    G <- add_vertices(graph = G, nv = nodes_add)
    n_a <- which(is.na(V(G)$name))
    l_n_a <- length(n_a)
    
    # Adding node times
    if(distribution == "continuous"){
      if(l_n_a > 0){
        # Their length of stay in the population
        V(G)$name[n_a] <- paste(i, "_", 1:l_n_a, sep = "")
        
        exit_time_new <- floor(rexp(l_n_a, mu)) + i + 1 # the "i" represents the years that have already passed
        names(exit_time_new) <- V(G)$name[n_a]
        exit_time <- c(exit_time, exit_time_new)
        
      }
      
      return(list(G, dynamic_edge_tracker, exit_time))
    }else{
      if(distribution == "discrete"){
        if(l_n_a > 0){
          # Their length of stay in the population
          V(G)$name[n_a] <- paste(i, "_", 1:l_n_a, sep = "")
          
          exit_time_new <- rep(Inf, l_n_a) # exit times will be decided discretely at each iteration
          names(exit_time_new) <- V(G)$name[n_a]
          exit_time <- c(exit_time, exit_time_new)
          
        }
        
        return(list(G, dynamic_edge_tracker, exit_time))
      }
    }
  }
  
  # Edge add
  Edge_Add <- function(G, n_partner, dynamic_edge_tracker, distribution, type){
    # G: Graph to be modified
    # n_partner: partners to potentially form steady relationship
    # dynamic_edge_tracker: Object tracking all edge information
    # distribution: Continuous or discrete
    # type: type of edge to add
    
    if(length(n_partner) > 1){
      
      if(length(n_partner)%%2 == 1){# if odd, add even number of edges
        steady_edge <- sample(x = n_partner,
                              size = length(n_partner) - 1,
                              replace = FALSE)
      }else{
        if(length(n_partner)%%2 == 0){# if even, add all partners
          steady_edge <- sample(x = n_partner,
                                size = length(n_partner),
                                replace = FALSE)
        }
      }
      
      G <- add_edges(graph = G,
                     edges = steady_edge)
      
      E_temp <- matrix(steady_edge, length(steady_edge)/2, 2, byrow = TRUE)
      
      if(type == "s"){
        if(distribution == "continuous"){
          temp <- cbind(E_temp, i, floor(rexp(dim(E_temp)[1], sigma)) + i + 1)
        }else{
          if(distribution == "discrete"){
            temp <- cbind(E_temp, i, Inf) # Separation is decided each iteration
          }
        }
      }else{
        if(type == "c"){
          temp <- cbind(E_temp, i, i + 1) # Separation is done on the next iteration
        }
      }
      
      dynamic_edge_tracker <- rbind(dynamic_edge_tracker, 
                                    cbind(temp, type))
    }
    
    return(list(G, dynamic_edge_tracker))
  }
  
  # Steady Edge Add
  Steady_Edge_Add <- function(G, dynamic_edge_tracker, i, rho, distribution){
    # G: Graph to be modified
    # dynamic_edge_tracker: Object tracking all edge information
    # i: Present iteration
    # rho: Probability of forming steady edge
    # distribution: Continuous or discrete
    
    if(distribution == "continuous"){
      steady_partner_rate <- rexp(length(V(G)), rho)
      names(steady_partner_rate) <- V(G)$name
      
      n_partner <- names(which(steady_partner_rate < 1))
      
    }else{
      if(distribution == "discrete"){
        steady_partner_rate <- rbinom(length(V(G)), 1, rho)
        names(steady_partner_rate) <- V(G)$name
        
        n_partner <- names(which(steady_partner_rate == 1))
      }
    }
    
    n_partner <- n_partner[n_partner %in% names(which(degree(G) == 0))]
    
    return(Edge_Add(G, n_partner, dynamic_edge_tracker, distribution, "s"))
  }
  
  # Casual Edge Add
  Casual_Edge_Add <- function(G, dynamic_edge_tracker, i, w_0, w_1, distribution){
    # G: Graph to be modified
    # dynamic_edge_tracker: Object tracking all edge information
    # i: Present iteration
    # w_0: Probability a single person forms a casual relationship
    # w_1: Probability a person in a relationship forms a casual relationship
    # distribution: Continuous or discrete
    
    single_nodes <- names(which(degree(G) == 0))
    steady_nodes <- names(which(degree(G) == 1))
    
    if(distribution == "continuous"){
      casual_single_rate <- rexp(length(single_nodes), w_0)
      names(casual_single_rate) <- single_nodes
      
      casual_steady_rate <- rexp(length(steady_nodes), w_1)
      names(casual_steady_rate) <- steady_nodes
      
      casual_rate <- sample(c(casual_single_rate, casual_steady_rate)) # randomly mixing all contacts
      
      n_partner <- names(which(casual_rate < 1))
    }else{
      if(distribution == "discrete"){
        casual_single_rate <- rbinom(length(single_nodes), 1, w_0)
        names(casual_single_rate) <- single_nodes
        
        casual_steady_rate <- rbinom(length(steady_nodes), 1, w_1)
        names(casual_steady_rate) <- steady_nodes
        
        casual_rate <- sample(c(casual_single_rate, casual_steady_rate)) # randomly mixing all contacts
        
        n_partner <- names(which(casual_rate == 1))
      }
    }
    
    return(Edge_Add(G, n_partner, dynamic_edge_tracker, distribution, "c"))
  }
  
  # Final Graph Puller
  Graph_Pull <- function(graph_tracker, dynamic_edge_tracker, g_union, lag, iter){
    # graph_tracker: tracks every iteration of the graph
    # dynamic_edge_tracker: Tracks all information on edges
    # g_union: number of iterations to union over
    # lag: Spacing between graphs
    # iter: Final iteration number in loop
    
    # Returning two time points based on g_union
    if(lag >= 0){
      
      G_union_list_0 <- graph_tracker[I(iter - g_union - g_union + 1):I(iter - g_union) - lag]
      G_union_list_1 <- graph_tracker[I(iter - g_union + 1):iter]
      
    }else{
      
      G_union_list_0 <- graph_tracker[I(iter - g_union + 1):iter]
      G_union_list_1 <- graph_tracker[I(iter - g_union + 1):iter]
    }
    
    G_union_0 <- G_union_list_0[[1]]
    G_union_1 <- G_union_list_1[[1]]
    
    if(g_union > 1){
      for(i in 2:g_union){
        G_union_0 <- igraph::union(G_union_0, G_union_list_0[[i]])
        G_union_1 <- igraph::union(G_union_1, G_union_list_1[[i]])
      }
    }
    
    G_union <- list(G_union_0, G_union_1, 
                    dynamic_edge_tracker) 
    
    return(G_union)
  }
  
  
  ## This section initializes graph and runs the mechanism
  G_0 <- make_empty_graph(n = num_nodes, directed = FALSE)
  V(G_0)$name <- paste(0, "_", 1:length(V(G_0)), sep = "") # allows me to delete particular nodes
  
  # initialing graph tracker
  graph_tracker <- list()
  # initializing edge tracker
  dynamic_edge_tracker <- matrix(NA, 0, 5)
  colnames(dynamic_edge_tracker) <- c("e1", "e2", "start", "end","type")
  # Each nodes continuous exit time
  exit_time <- ceiling(rexp(num_nodes, mu))
  names(exit_time) <- V(G_0)$name
  # Tracking number of nodes
  node_tracker <- NULL
  
  i = 1
  while(i <= iter | length(V(G_0)) != num_nodes){
    #print(c(i, length(V(G_0))))
    # Casual Edge Removal
    stage_0 <- Casual_Steady_Remove(G = G_0,
                                    dynamic_edge_tracker = dynamic_edge_tracker,
                                    i = i,
                                    sigma = sigma,
                                    distribution = distribution,
                                    type = "c")
    #print(0)
    
    
    # Steady Edge Removal
    stage_1 <- Casual_Steady_Remove(G = stage_0[[1]],
                                    dynamic_edge_tracker = stage_0[[2]],
                                    i = i,
                                    sigma = sigma, 
                                    distribution = distribution,
                                    type = "s")
    
    #print(1)
    # Node Removal
    stage_2 <- Node_Remove(G = stage_1[[1]],
                           dynamic_edge_tracker = stage_1[[2]],
                           exit_time = exit_time,
                           i = i,
                           mu = mu,
                           distribution = distribution)
    
    
    #print(2)
    # Add Nodes
    stage_3 <- Node_Add(G = stage_2[[1]],
                        dynamic_edge_tracker = stage_2[[2]],
                        exit_time = exit_time,
                        i = i,
                        mu = mu,
                        distribution = distribution)
    
    exit_time <- stage_3[[3]]
    
    #print(3)
    # Add Steady Edges
    stage_4 <- Steady_Edge_Add(G = stage_3[[1]],
                               dynamic_edge_tracker = stage_3[[2]],
                               i = i,
                               rho = rho,
                               distribution = distribution)
    
    #print(4)
    # Add Casual Edges
    stage_5 <- Casual_Edge_Add(G = stage_4[[1]],
                               dynamic_edge_tracker = stage_4[[2]],
                               i = i,
                               w_0 = w_0,
                               w_1 = w_1,
                               distribution = distribution)
    
    #print(5)
    # Recording Results
    G_0 <- graph_tracker[[i]] <- stage_5[[1]]
    dynamic_edge_tracker <- stage_5[[2]]
    node_tracker[i] <- length(V(G_0))
    
    i <- i + 1
    
  }
  
  # Only keeping relevant data from past 500 iterations
  thresh <- (i - 1) - 500
  dynamic_edge_tracker <- dynamic_edge_tracker[which(as.numeric(dynamic_edge_tracker[, "end"]) > thresh), ]
  iteration = i-1 # Final iteration pulled from while loop
  dynamic_edge_tracker <- cbind(dynamic_edge_tracker, iteration)
  
  return(Graph_Pull(graph_tracker = graph_tracker,
                    dynamic_edge_tracker = dynamic_edge_tracker,
                    g_union = g_union,
                    lag = lag,
                    iter = iteration)) # we stop before ith iteration
}
 
