setwd("~/Documents/PhD/Chapter 4/FINAL DATA AND CODE/RNA Seq")

library(DESeq2)
library(readr)
library(dplyr)
library(tidyr)
library(tibble)

# Get transcripts with total TPM over 22
# ChrisV used TPM sum > 30 over 79 samples = 0.38 TPM/sample
# We have 58 samples. * 0.38TPM = 22
tpms <- read_csv("1-Process-HTSeq/output/TPM.csv") %>%
    column_to_rownames("transcript_id") %>%
    filter(rowSums(.) > 22) %>%
    rownames()

# Get information about each transcript ID
gene_to_go <- read_tsv("3-Reference_to_GO/output/gene_to_go.tsv")

## DEGs

# run dds_data.r which we used to create a DESeq dataset from our data
source("2-Differential_gene_expression/dds_data.r")

# levels(metadata$Group)
# Add columns indicating whether the sample is in or not in each group
metadata  <- metadata %>% mutate(
    is_mangrove_reef = ifelse(Group == "wild_mangrove", 'excluded', Group == "mangrove_reef"),
    is_reef_reef = ifelse(Group == "wild_reef", 'excluded', Group == "reef_reef"),
    is_wild_mangrove = ifelse(Group == "mangrove_reef", 'excluded', Group == "wild_mangrove"),
    is_wild_reef = ifelse(Group == "reef_reef", 'excluded', Group == "wild_reef")
) %>% mutate(
    is_mangrove_reef = as.factor(as.character(is_mangrove_reef)),
    is_reef_reef = as.factor(as.character(is_reef_reef)),
    is_wild_mangrove = as.factor(as.character(is_wild_mangrove)),
    is_wild_reef = as.factor(as.character(is_wild_reef))
)

# use new metadata in dds
dds <- DESeqDataSetFromMatrix(
  countData = countData,
  colData = metadata, design = ~Temperature
)

dds <- DESeq(dds)

# Transcripts upregulated at 33 (DESeq2 Results)
up33v30 <- results(dds, contrast = c("Temperature", "33", "30")) %>%
    as.data.frame() %>%
    filter(log2FoldChange > 1.99, pvalue < 0.005) %>%
    rownames()

# Transcripts upregulated at 36 (DESeq2 Results)
up36v30 <- results(dds, contrast = c("Temperature", "36", "30")) %>%
    as.data.frame() %>%
    filter(log2FoldChange > 1.99, pvalue < 0.005) %>%
    rownames()

# Transcripts upregulated at 33 or 36
up33or36 <- union(up33v30, up36v30)

## Function to find front-loaded genes

# target = which one might do we want to check for front loading
# other = front loaded compared to what?
find_front_loaded <- function(target) {
    message("Checking for front loading in ", target, " compared to all others")

    design(dds) <- as.formula(paste0("~is_", target))
    dds <- DESeq(dds)

    # Transcripts upregulated in `target` vs all others at baseline (30º) (DESeq2 Results)
    up_target <- results(dds, contrast = c(paste0("is_", target), 'TRUE', 'FALSE')) %>%
        as.data.frame() %>%
        filter(log2FoldChange > 1.99, pvalue < 0.005) %>%
        rownames()

    # Transcripts upregulated at (33 or 36) and in `target` vs all others
    up <- intersect(up33or36, up_target)

    ## TPMs

    # Transcripts with total TPM over threshold and upregulated at (33 or 36) and upregulated in all others vs `target`
    front_loaded <- intersect(up, tpms)

    message("Found ", length(front_loaded), " front loaded transcripts")

    gene_to_go %>%
        filter(gene_id %in% front_loaded) %>%
        write.csv(paste0("2-Differential_gene_expression/output/front_loading/vs_all_2_", target, ".csv"))
}

groups <- c(
    "reef_reef",
    "wild_mangrove",
    "wild_reef",
    "mangrove_reef"
)
for (group in groups) {
    find_front_loaded(group)
}
