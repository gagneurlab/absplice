library(dplyr)
library(FRASER)
library(data.table)

saveFdsAsCountTable <- function(fds, min_read) {
  # junction info
  junction_dt <- as.data.table(granges(rowRanges(fds)))
  junction_dt$width <- NULL # remove uncessary columns
  colnames(junction_dt) <- c("Chromosome", "Start", "End", "Strand") # rename to fulfill the format https://gitlab.cmm.in.tum.de/gagneurlab/count_table
  junction_dt$Start <- as.integer(junction_dt$Start - 1) # change 1 based to 0 based. see https://www.biostars.org/p/84686/
  print(table(junction_dt$Strand))
  junction_dt$Strand <- gsub("[*]",".", junction_dt$Strand) # replace * with . 
  print(table(junction_dt$Strand))
  
  # junction count
  count_dt <- as.matrix(counts(fds))
  count_csv <- cbind(junction_dt, count_dt)
  
  # remove sparse rows
  count_csv <- count_csv[rowSums(count_dt) >= min_read, ]
  return(count_csv)  
}

working_dir <- snakemake@params[['fraser_working_dir']]
analysis_name <- paste0(snakemake@params[['analysis']])

fds <- loadFraserDataSet(dir = working_dir, name = analysis_name)
count_csv <- saveFdsAsCountTable(fds, min_read = 0)
fwrite(count_csv, file = snakemake@output[['count_table']])