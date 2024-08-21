# Converts HTSeq output files to CSV for use in analysis


library(readr)
library(dplyr)
library(tidyr)

setwd("~/Documents/PhD/Chapter 4/Data/R analysis/RNA-seq/P. acuta new genome")

metadata <- read_csv("metadata.csv") %>%
    # Add unique identifier column
    mutate(Identifier = paste0("t", Tank, Group, SampleID))

data <- NULL
diagnostics <- NULL
for (i in seq_len(nrow(metadata))) {
    identifier <- metadata[[i, "Identifier"]]

    filepath <- paste0("1-Process-HTSeq/HTSeq_data/", metadata[i, "FileName"], ".htseq.out")
    message(filepath)
    file <- read_tsv(filepath, col_names = c("Gene", identifier), show_col_types = FALSE)

    if(nrow(file)==0){
        message("!!!File is empty!!!")
        next
    }

    if (i == 1) {
            # First run, put the data in `data` and `diagnostics`
        data <- file %>% filter(!startsWith(Gene, "__"))
        diagnostics <- file %>% filter(startsWith(Gene, "__"))
    } else {
        # Second and subsequent runs, join with the existing data
        data <-
            left_join(
                data,
                file %>% filter(!startsWith(Gene, "__")),
                "Gene"
            )
        diagnostics <-
            left_join(
                diagnostics,
                file %>% filter(startsWith(Gene, "__")),
                "Gene"
            )
    }
}

# Find NAs
# which(is.na(data), arr.ind = TRUE)

write.csv(data, "1-Process-HTSeq/out/HTSeqCounts.csv", row.names = FALSE)
write.csv(diagnostics, "1-Process-HTSeq/out/HTSeqDiagnostics.csv", row.names = FALSE)
