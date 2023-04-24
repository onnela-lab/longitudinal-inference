####################################################################################################################################
# This script calculates the edge index from a Network when given a sequence of edges as its input. If the edge sequence
# is one node, then all indicies involving that edge are returned. 
####################################################################################################################################


edge_index <- function(G, 
                       edge_sequence,
                       edge_list = NULL){
  # G:
    # network graph
  # edge_sequence:
    # sequence of node pairs that captures the edges
    # c(1,2,5,3,8,7) represents edges 1-2, 5-3, 8-7
    # or singular node that identifies all edges said node is a part of
  # edge_list:
    # An edge list can be inputted instead of a graph
    # if needed.
  
  if(is.null(edge_list)){
    full_edge_list <- get.edgelist(G)
  }else{
    full_edge_list <- edge_list
  }
  saved_index <- NULL
  
  if(length(edge_sequence) == 1){
    edge <- edge_sequence
    n <- which((full_edge_list[, 1] == edge[1])|
                 (full_edge_list[, 2] == edge[1]))
    saved_index <- n
    
  }else{
    if(length(edge_sequence) > 1){
      for(i in seq(2, length(edge_sequence), by = 2)){
        edge <- edge_sequence[I(i-1):i]
        n <- which((full_edge_list[, 1] == edge[1] & full_edge_list[, 2] == edge[2] )|
                     (full_edge_list[, 2] == edge[1] & full_edge_list[, 1] == edge[2] ))
        
        # If there are multi edges check if there is a third
        # vector on the matrix, if there is a third vector
        # select the edge that has a 0 next to it. 
        # if no third vector or multiple 0s, select at random.
        if(!is.null(edge_list)){
          if(dim(edge_list)[2] > 2 & length(n) > 1){
            n <- n[which(edge_list[n, 3] == 0)]
          }
        }
        
        if(length(n) > 1){ # sample n such that we didn't already sample it
          n <- sample(n[which(!n %in% saved_index)], 1)
        }else{ #sample is weird when inputting one number to sample from
          n <- n
        }
        
        if(length(n) == 1){
          saved_index[i/2] <- n
        }else{
          saved_index[i/2] <- NA
        }
      }
    }
  }
  
  return(saved_index)
}
