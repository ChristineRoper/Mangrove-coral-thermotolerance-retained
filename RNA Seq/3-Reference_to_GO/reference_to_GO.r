
# Unload and reload dplyr. Because the global environment is shared in R, and
# libraries can overwrite things, something keeps replacing dplyr functions with
# versions that do differnt things, so the code breaks. We need to force them
# back to the functions we expect.
unloadNamespace("tidyr")
unloadNamespace("dplyr")
library(dplyr)
library(tidyr)
library(readr)
library(stringr)

## Objective: Use the reference genome files to map genes -> transcript ID -> GO name


setwd("~/Documents/PhD/Chapter 4/FINAL DATA AND CODE/RNA Seq")

## This file tells us which GO names go with which transcript IDs (both in reference genome)
## file name: Pocillopora_acuta_HIv2.genes.EggNog_results.txt

transcript_to_go <- read_tsv("reference genome/Pocillopora_acuta_HIv2.genes.EggNog_results.txt") %>%
    # rename columns
    mutate(transcript_id = `#query`) %>% # rename "#query" column to "transcript_ID"
    # Only keep these four columns
    select(transcript_id, GOs, Description, Preferred_name)

## This file tells us which transcript_id goes with which gene_id
## Pocillopora_acuta_HIv2.genes.gtf

gene_to_transcript <- read_tsv(
    "reference genome/Pocillopora_acuta_HIv2.genes.gtf",
    skip = 1, # skip row one, contains no data
    col_names = c( #columns do not have names so here we are giving them names
        "SequenceID", "Source",
        "FeatureType", "FeatureStart", "FeatureEnd",
        "Score", "Stand", "Phase",
        "Attributes"
    )
) %>%
    # Only keep the "transcripts" because they have the transcript_id and gene_id
    filter(FeatureType == "transcript") %>%
    # We only need info in the "Attributes" column
    select(Attributes) %>%
    # Split it to get the gene_id and transcript_id, if there are more attributes, drop them. 
    separate_wider_delim(Attributes, delim = regex("; | "), names = c(NA, "gene_id", NA, "transcript_id"), too_many = "drop") %>%
    # Strip the quotemarks "Pocillopora_acuta_HIv2___TS.g24026.t1" => Pocillopora_acuta_HIv2___TS.g24026.t1
    mutate(gene_id = gsub('"', "", gene_id), transcript_id = gsub('"', "", transcript_id))


### Joining gene_to_transcript with transcript_to_go we get gene_id_to_go

left_join(gene_to_transcript, transcript_to_go, by = "transcript_id") %>%
    select(gene_id, transcript_id, GOs, Description, Preferred_name) %>%
    write_tsv("3-Reference_to_GO/output/gene_to_go.tsv")

### Do it again in the format that TopGO likes

left_join(gene_to_transcript, transcript_to_go, by = "transcript_id") %>%
    filter(!is.na(GOs))%>% # only keep the genes with GO terms
    filter(GOs != '-')%>% # don't keep any that are dashed (no terms)
    select(gene_id, GOs) %>%
    write_tsv("3-Reference_to_GO/output/gene_to_TopGO.tsv", col_names = FALSE)
