setwd("RNA Seq")

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

design(dds) <- ~TempGroup
dds <- DESeq(dds)

## Function to find front-loaded genes

# target = which one might do we want to check for front loading
# other = front loaded compared to what?
find_front_loaded <- function(target, other) {
    message("Checking for front loading in ", target, " compared to ", other)

    # Transcripts upregulated at 33 (DESeq2 Results)
    up33v30 <- results(dds, contrast = c("TempGroup", paste0("33_", other), paste0("30_", other))) %>%
        as.data.frame() %>%
        filter(log2FoldChange > 1.99, pvalue < 0.005) %>%
        rownames()

    # Transcripts upregulated at 36 (DESeq2 Results)
    up36v30 <- results(dds, contrast = c("TempGroup", paste0("36_", other), paste0("30_", other))) %>%
        as.data.frame() %>%
        filter(log2FoldChange > 1.99, pvalue < 0.005) %>%
        rownames()

    # Transcripts upregulated at 33 or 36
    up33or36 <- union(up33v30, up36v30)

    # Transcripts upregulated in `target` vs `other` at baseline (30º) (DESeq2 Results)
    up_target <- results(dds, contrast = c("TempGroup", paste0("30_", target), paste0("30_", other))) %>%
        as.data.frame() %>%
        filter(log2FoldChange > 1.99, pvalue < 0.005) %>%
        rownames()

    # Transcripts upregulated at (33 or 36) and in `target` vs `other`
    up <- intersect(up33or36, up_target)

    ## TPMs

    # Transcripts with total TPM over threshold and upregulated at (33 or 36) and upregulated in `other` vs `target`
    front_loaded <- intersect(up, tpms)

    message("Found ", length(front_loaded), " front loaded transcripts")

    gene_to_go %>%
        filter(gene_id %in% front_loaded) %>%
        write.csv(paste0("2-Differential_gene_expression/output/front_loading/", target, "_vs_", other, ".csv"))
}

resultsNames(dds)
groups <- c(
    "reef_reef",
    "wild_mangrove",
    "wild_reef",
    "mangrove_reef"
)
for (groupA in groups) {
    for (groupB in groups) {
        if (groupA == groupB) next
        find_front_loaded(groupA, groupB)
    }
}
