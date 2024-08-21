library(DESeq2)
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(tibble)

setwd("~/Documents/PhD/Chapter 4/Data/R analysis/RNA-seq/P. acuta new genome")

############### Overview of gene expression using DESeq2 ###################

metadata <- read_csv("metadata.csv") %>%
  # Add Origin and Destination
  separate(Group, c("Origin", "Destination"), remove = FALSE) %>%
  # Add unique identifier column
  mutate(Identifier = paste0("t", Tank, Group, SampleID)) %>%
  # Wild ones came from their destination
  mutate(Moved = Origin != "wild", Origin = ifelse(Origin == "wild", Destination, Origin)) %>%
  # Make a TempGroup column
  tidyr::unite(TempGroup, Temperature, Group, remove = FALSE) %>%
  # make factors
  mutate(
    Identifier = factor(Identifier),
    Tank = factor(Tank),
    Temperature = factor(Temperature),
    Group = factor(Group),
    TempGroup = factor(TempGroup),
    SampleID = factor(SampleID),
    Origin = factor(Origin),
    Destination = factor(Destination)
  )

# HTSeq
#######

# input HTSeq gene count matrix
countData <- as.matrix(read.csv("1-Process-HTSeq/out/HTSeqCounts.csv", row.names = "Gene"))

# check names match in files
all(metadata$Identifier %in% colnames(countData))
all(metadata$Identifier == colnames(countData))

# # Check for levels with no samples
# m1 <- model.matrix(~ Group * Temperature, metadata)
# colnames(m1)
# unname(m1)
# # Should all be false
# apply(m1, 2, function(x) all(x==0))

dds <- DESeqDataSetFromMatrix(
  countData = countData,
  colData = metadata, design = ~TempGroup
)
