# install.packages("BiocManager")
# BiocManager::install("GSEABase")
# BiocManager::install("GO.db")

library(readr)
library(dplyr)
library(tibble)
library(GSEABase)
library(ggplot2)

slim <- getOBOCollection("5-Goslim/goslim_generic.obo")
slim_ids <- GSEABase::ids(slim)

# The oposite of goSlim - given slim, find genes
goFat <- function(slim, genes, GO_ontology) {
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
    "36_mangrove_reef_vs_reef_reef_p0.05.txt",
    "36_mangrove_reef_vs_wild_mangrove_p0.05.txt",
    "36_mangrove_reef_vs_wild_reef_p0.05.txt",
    "36_reef_reef_vs_wild_mangrove_p0.05.txt",
    "36_reef_reef_vs_wild_reef_p0.05.txt",
    "36_wild_mangrove_vs_wild_reef_p0.05.txt",
    "mangrove_reef_30_vs_36_p0.05.txt",
    "reef_reef_30_vs_36_p0.05.txt",
    "wild_mangrove_30_vs_36_p0.05.txt",
    "wild_reef_30_vs_36_p0.05.txt",
    # Added 2024-08-21
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
    functions <- read_tsv(paste0("4-Functional_annotation/out/combined_", direction, ".tsv"))


    # output <- NULL
    for (current_file in files) {
        message(current_file)

        GO_IDs <- filter(functions, file == current_file, regulation == direction)$GO.ID
        print(GO_IDs)
        collection <- GOCollection(GO_IDs)

        for (ontology in c("MF", "BP", "CC")) {
            message(ontology)

            results <- goSlim(collection, slim, ontology) %>%
                rownames_to_column("GO") %>%
                mutate(
                    File = current_file,
                    direction = direction,
                    Ontology = ontology
                ) %>%
                rowwise() %>%
                mutate(genes = paste(goFat(GO, GO_IDs, ontology), collapse = ", "))
            # output <- rbind(output, results)
            results %>%
                filter(Count > 0) %>%
                write_tsv(paste0("5-Goslim/out/", direction, ".tsv"), append = TRUE)
        }
    }

}
