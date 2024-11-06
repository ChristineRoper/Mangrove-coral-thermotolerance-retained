# BiocManager::install("topGO")

library("topGO")
library("data.table")
library("ggplot2")
library("RColorBrewer")
library("ggVennDiagram")
library("dplyr")
library(stringr)
library(readr)

setwd("~/Documents/PhD/Chapter 4/FINAL DATA AND CODE/RNA Seq")

############### GO enrichment using topGO ###################

geneID2GO <- readMappings(file = "3-Reference_to_GO/output/gene_to_TopGO.tsv")
# all possible gene names with GO annotation
geneUniverse <- names(geneID2GO)
message("Number of genes with GO annotation: ", length(geneUniverse))

files <- c(
  "by-group/wild_mangrove_vs_wild_reef_p0.05.txt",
  "by-group/reef_reef_vs_wild_reef_p0.05.txt",
  "by-group/reef_reef_vs_wild_mangrove_p0.05.txt",
  "by-group/mangrove_reef_vs_wild_reef_p0.05.txt",
  "by-group/mangrove_reef_vs_wild_mangrove_p0.05.txt",
  "by-group/mangrove_reef_vs_reef_reef_p0.05.txt",
  "by-temperature/36_vs_33_p0.05.txt",
  "by-temperature/36_vs_30_p0.05.txt",
  "by-temperature/33_vs_30_p0.05.txt",
  "group-by-temperature/30_mangrove_reef_vs_reef_reef_p0.05.txt",
  "group-by-temperature/30_mangrove_reef_vs_wild_mangrove_p0.05.txt",
  "group-by-temperature/30_mangrove_reef_vs_wild_reef_p0.05.txt",
  "group-by-temperature/30_reef_reef_vs_wild_mangrove_p0.05.txt",
  "group-by-temperature/30_reef_reef_vs_wild_reef_p0.05.txt",
  "group-by-temperature/30_wild_mangrove_vs_wild_reef_p0.05.txt",
  "group-by-temperature/33_mangrove_reef_vs_reef_reef_p0.05.txt",
  "group-by-temperature/33_mangrove_reef_vs_wild_mangrove_p0.05.txt",
  "group-by-temperature/33_mangrove_reef_vs_wild_reef_p0.05.txt",
  "group-by-temperature/33_reef_reef_vs_wild_mangrove_p0.05.txt",
  "group-by-temperature/33_reef_reef_vs_wild_reef_p0.05.txt",
  "group-by-temperature/33_wild_mangrove_vs_wild_reef_p0.05.txt",
  "group-by-temperature/36_mangrove_reef_vs_reef_reef_p0.05.txt",
  "group-by-temperature/36_mangrove_reef_vs_wild_mangrove_p0.05.txt",
  "group-by-temperature/36_mangrove_reef_vs_wild_reef_p0.05.txt",
  "group-by-temperature/36_reef_reef_vs_wild_mangrove_p0.05.txt",
  "group-by-temperature/36_reef_reef_vs_wild_reef_p0.05.txt",
  "group-by-temperature/36_wild_mangrove_vs_wild_reef_p0.05.txt"
)


for (direction in c("up", "down")) {
  for (file in files) {
    message(direction, " ", file)

    filepath <- paste("2-Differential_gene_expression/output/pairwise", file, sep = "/")
    DE <- read.delim(filepath, sep = "\t", header = TRUE) # read differentially expressed genes

    if (direction == "up") {
      DE <- filter(DE, stat > 0) # filter to only include upregulated genes (ie. if the value is > 0)
    }
    if (direction == "down") {
      DE <- filter(DE, stat < 0) # filter to only include downregulated genes (ie. if the value is < 0)
    }

    DE$gene_names <- row.names(DE) # make a column called "gene_names" from the row name
    DE_genes <- DE$gene_names
    message("Number of genes in file: ", length(DE_genes))

    # keep the differentially expressed genes which are in geneUniverse
    keep <- DE_genes %in% geneUniverse
    keep <- which(keep == TRUE)
    DE_genes <- DE_genes[keep]
    message("Number of genes in file with GO annotation: ", length(DE_genes))
    # TODO: Look at the ones without annotation and try to find out why

    # make named vector list of factors showing which are genes of interest
    geneList <- factor(as.integer(geneUniverse %in% DE_genes)) # creates list of all genes represented by 1L (not in the DEGs list) or 2L (in the DEGs list) which is required for TopGo
    names(geneList) <- geneUniverse

    options(scipen = 999) # change format to non-scientific digits

    # BP
    # create topgo data object
    myGOdata <- new("topGOdata", description = "My project", ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
    # test for significance
    # run weighted algorithm as classic doesnt take into consideration GO hierarchy so could overrepresent enrichment
    resultFisher <- runTest(myGOdata, algorithm = "weight01", statistic = "fisher")
    # generate a table of results using Gentable function
    allGO <- usedGO(object = myGOdata)
    BP <- GenTable(myGOdata, weightFisher = resultFisher, orderBy = "resultsFisher", ranksOf = "weightFisher", topNodes = length(allGO))
    BP$ontology <- "BP"

    # MF
    myGOdata <- new("topGOdata", description = "My project", ontology = "MF", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
    resultFisher <- runTest(myGOdata, algorithm = "weight01", statistic = "fisher")
    allGO <- usedGO(object = myGOdata)
    MF <- GenTable(myGOdata, weightFisher = resultFisher, orderBy = "resultFisher", ranksOf = "weightFisher", topNodes = length(allGO))
    MF$ontology <- "MF"

    # CC
    myGOdata <- new("topGOdata", description = "My project", ontology = "CC", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
    resultFisher <- runTest(myGOdata, algorithm = "weight01", statistic = "fisher")
    allGO <- usedGO(object = myGOdata)
    CC <- GenTable(myGOdata, weightFisher = resultFisher, orderBy = "resultFisher", ranksOf = "weightFisher", topNodes = length(allGO))
    CC$ontology <- "CC"

    # Creates table which combines weightFisher values from BP, CC and MF within our cut-off (< 0.001)
    out <- rbind(BP, CC, MF)
    out_significant <- subset(out, out$weightFisher < 0.001)
    
    write.table(out_significant, paste0("4-Functional_annotation/output/", direction, "-", file), row.names = FALSE, sep = "\t", quote = FALSE)
  }
}

# combine the above outputs for each direction into one file
for (direction in c("up", "down")) {
  combined <- NULL
  for (file in files) {
    file_data <- read_tsv(paste0("4-Functional_annotation/output/", direction, "-", file))
    file_data$file <- file
    file_data$regulation <- direction

    combined <- rbind(combined, file_data)
  }
  write_tsv(combined, paste0("4-Functional_annotation/output/combined_", direction, ".tsv"))
}
 