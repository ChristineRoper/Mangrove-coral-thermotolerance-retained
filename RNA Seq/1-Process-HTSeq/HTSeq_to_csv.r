# Combines HTSeq output files and converts to csv for use in analysis

library(readr)
library(dplyr)
library(tidyr)

# input: metadata.csv and raw HTSeq_data 
# output: HTSeqCounts.csv and HTSeqDiagnostics.csv

setwd("~/Documents/PhD/Chapter 4/FINAL DATA AND CODE/RNA Seq")

# Read metadata
metadata <- read_csv("metadata.csv") %>%
    # Add unique identifier column by combining tank, group and sample ID
    mutate(Identifier = paste0("t", Tank, Group, SampleID))

data <- NULL
diagnostics <- NULL # move diagnostics data from HTSeq output files and combine in separate csv

for (i in seq_len(nrow(metadata))) {
    identifier <- metadata[[i, "Identifier"]] # get identifier for each row

    filepath <- paste0("1-Process-HTSeq/HTSeq_data/", metadata[i, "FileName"], ".htseq.out") # construct the file path
    message(filepath)
    file <- read_tsv(filepath, col_names = c("Gene", identifier), show_col_types = FALSE) # read the HTSeq file

    # if the file is empty, move to the next one (for any samples that failed QC)
    if (nrow(file) == 0) {
        message("!!!File is empty!!!")
        next
    }

    file_data <- file %>% filter(!startsWith(Gene, "__")) # select all data that doesn't start with "__" (this is the diagnostic data)
    file_diagnostics <- file %>% filter(startsWith(Gene, "__")) # select the diagnostic data (starts with "__")

    if (i == 1) {
        # First file/run put the data in 'data' and 'diagnostics'
        data <- file_data
        diagnostics <- file_diagnostics
    } else {
        # Second and subsequent runs, join with the existing data from the previous files
        data <- left_join(data, file_data, "Gene") 
        diagnostics <- left_join(diagnostics, file_diagnostics, "Gene")
    }
}


write.csv(data, "1-Process-HTSeq/output/HTSeqCounts.csv", row.names = FALSE)
write.csv(diagnostics, "1-Process-HTSeq/output/HTSeqDiagnostics.csv", row.names = FALSE)
