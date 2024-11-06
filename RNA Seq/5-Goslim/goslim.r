# install.packages("BiocManager")
# BiocManager::install("GSEABase")
# BiocManager::install("GO.db")

library(readr)
library(dplyr)
library(tibble)
library(GSEABase)
library(ggplot2)

setwd("~/Documents/PhD/Chapter 4/FINAL DATA AND CODE/RNA Seq")

# Load the GO Slim file
slim <- getOBOCollection("5-Goslim/goslim_generic.obo")
slim_ids <- GSEABase::ids(slim) # gets slim IDs

# The oposite of GOSlim - given slim, this function finds the genes
# Based on the GOSlim function, reverse engineer the process to give us the genes
goWide <- function(slim, genes, GO_ontology) {
    GOTERM <- getAnnMap("TERM", "GO")
    terms <- mget(slim_ids, GOTERM, ifnotfound = NA)
    terms <- terms[vapply(terms, Ontology, character(1)) == GO_ontology]
    slim_ids_o <- names(terms)

    ## Use GO_ontology to find the required offspring
    OFFSPRING <- switch(GO_ontology,
        MF = getAnnMap("MFOFFSPRING", "GO"),
        BP = getAnnMap("BPOFFSPRING", "GO"),
        CC = getAnnMap("CCOFFSPRING", "GO"),
        stop("GO_ontology must be 'MF', 'BP', or 'CC'")
    )

    slim_to_gene <- mget(slim_ids_o, OFFSPRING, ifnotfound = NA)
    slim_to_gene <- slim_to_gene[!is.na(slim_to_gene)]
    gene_to_slim <- revmap(slim_to_gene)[genes]
    gene_to_slim <- gene_to_slim[!vapply(gene_to_slim, is.null, logical(1))]

    if (length(gene_to_slim) == 0) {
        return(c())
    }

    return(names(gene_to_slim)[sapply(gene_to_slim, function(x) any(grepl(slim, unlist(x))))])
}

files <- c(
    "group-by-temperature/36_mangrove_reef_vs_reef_reef_p0.05.txt",
    "group-by-temperature/36_mangrove_reef_vs_wild_mangrove_p0.05.txt",
    "group-by-temperature/36_mangrove_reef_vs_wild_reef_p0.05.txt",
    "group-by-temperature/36_reef_reef_vs_wild_mangrove_p0.05.txt",
    "group-by-temperature/36_reef_reef_vs_wild_reef_p0.05.txt",
    "group-by-temperature/36_wild_mangrove_vs_wild_reef_p0.05.txt",
    "group-by-temperature/mangrove_reef_30_vs_36_p0.05.txt",
    "group-by-temperature/reef_reef_30_vs_36_p0.05.txt",
    "group-by-temperature/wild_mangrove_30_vs_36_p0.05.txt",
    "group-by-temperature/wild_reef_30_vs_36_p0.05.txt",
    "by-group/wild_mangrove_vs_wild_reef_p0.05.txt",
    "by-group/reef_reef_vs_wild_reef_p0.05.txt",
    "by-group/reef_reef_vs_wild_mangrove_p0.05.txt",
    "by-group/mangrove_reef_vs_wild_reef_p0.05.txt",
    "by-group/mangrove_reef_vs_wild_mangrove_p0.05.txt",
    "by-group/mangrove_reef_vs_reef_reef_p0.05.txt",
    "by-temperature/36_vs_33_p0.05.txt",
    "by-temperature/36_vs_30_p0.05.txt",
    "by-temperature/33_vs_30_p0.05.txt"
)

for (direction in c("up", "down")) {
    functional_annotations <- read_tsv(paste0("4-Functional_annotation/output/combined_", direction, ".tsv")) # read summary tsv created in the previous step (either up or down direction)

    for (current_file in files) {
        message(direction, ' ', current_file) # tells you the file name

        GO_IDs <- filter(functional_annotations, file == current_file, regulation == direction)$GO.ID # filtering the combined functional annotation files to get the GO ID's for only the comparisons we're interested in
        collection <- GOCollection(GO_IDs) # making a collection required for GOSlim

        for (ontology in c("MF", "BP", "CC")) {
            message(ontology)

            # run GoSlim
            results <- goSlim(collection, slim, ontology) %>%
                rownames_to_column("GO") %>%
                # add these three columns
                mutate(
                    File = current_file,
                    direction = direction,
                    Ontology = ontology
                ) %>%
                rowwise() %>%
                mutate(GOTerms = paste(goWide(GO, GO_IDs, ontology), collapse = ", ")) %>% # adding new column for the GO terms
                filter(Count > 0) %>%
                write_tsv(paste0("5-Goslim/output/", direction, ".tsv"), append = TRUE)
        }
    }
}
