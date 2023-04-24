######################################################################
# This script contains all string processing functions
######################################################################
library(stringr)



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
  
  for(j in 1:length(iter_removed_1)){ # removing potential "_" in beginning of string
    if(substr(iter_removed_1[j], 1, 1) == "_"){
      iter_removed_1[j] <- substr(iter_removed_1[j],
                                  2,
                                  nchar(iter_removed_1[j]))
    }
  }
  
  # Just searching for edgelist lists
  iter_removed_edge <- unique(grep("edgeList", iter_removed_1, value = TRUE))
  
  # returning 1 iteration of all of the parameters independent of time varying
  return(iter_removed_edge)
  
}

# This function obtains the value parameter from the string
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

# This function obtains the numbers from our joint distribution strings
joint_processor <- function(str,
                            param){
  
  # str:
  # String we are processing 
  # param:
  # parameter we are matching on
  
  # Obtaining whole numbers first (discrete distribution)
  if(param %in% c("k_max", "g_union", "node", "lag", "size", 
                  "std")){# change std to decimal at some point
    i_results <- str_extract(str,
                             paste(param,
                                   ".{1}\\d*_", sep = "_"))
  }
  
  
  # Obtaining numbers with decimals next (chosen from continuous
  # distribution so "." guaranteed)
  if(param %in% c("gamma", "p_0", "p_cas", "migration", "m",
                  "mu", "rho", "sigma",  "w_0", "w_1", "theta",
                  "p_r", "p_d", "delta")){
    i_results <- str_extract(str,
                             paste(param,
                                   "\\d*.{0,1}\\d*_", sep = "_"))
  }
  
  # Removing "_" from end of string for processing later
  return(substr(i_results,
                1,
                nchar(i_results)-1))
  
}
