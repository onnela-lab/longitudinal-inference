################################################################################################
# This script takes in two nX2 matricies and counts how many rows in mat_2 are currently in 
# mat_1 and returns those edges. In essence, how many edges are maintained from 1 ---> 2
################################################################################################

common_edges <- function(mat_1, mat_2){
  # mat_1
    # nx2 matrix
  # mat_2
    # nx2 matrix
  
  mat_1_0 <- paste(mat_1[, 1], ".", mat_1[, 2], sep = "")
  mat_1_1 <- paste(mat_1[, 2], ".", mat_1[, 1], sep = "")
  
  mat_2_01 <- paste(mat_2[, 1], ".", mat_2[, 2], sep = "")
  
  if(dim(mat_1)[1] == 0 | dim(mat_2)[1] == 0){
    n_0 <- NULL
  }else{
    n_0 <- which(mat_2_01 %in% mat_1_0 | mat_2_01 %in% mat_1_1)
  }
  
  if(length(n_0) == 0){
    c_edges <- matrix(".", 0, 2)
  }else{
    c_edges <- mat_2[n_0, ,drop = FALSE]
  }
  
  return(c_edges)
}

