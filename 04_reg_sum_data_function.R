########################################################################
# This script combines all raw files into the reference table
########################################################################

# Parse in command line output:
cmdArgs <- commandArgs(trailingOnly = TRUE) 
z <- as.character(cmdArgs[1])

wd <- file.path("/home/ot25/Research/JP/HIV/Results/Reference_Table",
                z,
                "Raw_Files")

setwd(wd)

all_files <- list.files()


file <- paste(1, 
              "_",
              z,
              "_ABC_stat_table.csv",
              sep = "")

final <- read.csv(file)

for(i in 2:length(all_files)){
  print(i/length(all_files))
  file <- paste(i, 
                "_",
                z,
                "_ABC_stat_table.csv",
                sep = "")
  
  final <- rbind(final, read.csv(file))
  

}

wd <- file.path("/home/ot25/Research/JP/HIV/Results/Reference_Table",
                z)

setwd(wd)

file_name <- paste(z, "ABC_stat_table.csv", sep = "_")
write.csv(final, 
          file = file_name, 
          row.names = FALSE)
