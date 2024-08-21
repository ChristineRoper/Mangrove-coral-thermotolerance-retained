# BiocManager::install("topGO")

library("topGO")
library("data.table")
library("ggplot2")
library("RColorBrewer")
library("ggVennDiagram")
library("dplyr")
library(stringr)
library(readr)

############### GO enrichment using topGO ###################

geneID2GO <- readMappings(file = "3-Reference_to_GO/out/gene_to_go.tsv")
# all possible gene names with GO annotation
geneUniverse <- names(geneID2GO)
message("Number of genes with GO annotation: ", length(geneUniverse))

files <- c(
  "30_mangrove_reef_vs_reef_reef_p0.05.txt",
  "30_mangrove_reef_vs_wild_mangrove_p0.05.txt",
  "30_mangrove_reef_vs_wild_reef_p0.05.txt",
  "30_reef_reef_vs_wild_mangrove_p0.05.txt",
  "30_reef_reef_vs_wild_reef_p0.05.txt",
  "30_wild_mangrove_vs_wild_reef_p0.05.txt",
  "33_mangrove_reef_vs_reef_reef_p0.05.txt",
  "33_mangrove_reef_vs_wild_mangrove_p0.05.txt",
  "33_mangrove_reef_vs_wild_reef_p0.05.txt",
  "33_reef_reef_vs_wild_mangrove_p0.05.txt",
  "33_reef_reef_vs_wild_reef_p0.05.txt",
  "33_wild_mangrove_vs_wild_reef_p0.05.txt",
  "36_mangrove_reef_vs_reef_reef_p0.05.txt",
  "36_mangrove_reef_vs_wild_mangrove_p0.05.txt",
  "36_mangrove_reef_vs_wild_reef_p0.05.txt",
  "36_reef_reef_vs_wild_mangrove_p0.05.txt",
  "36_reef_reef_vs_wild_reef_p0.05.txt",
  "36_wild_mangrove_vs_wild_reef_p0.05.txt",
  "mangrove_reef_33_vs_30_p0.05.txt",
  "mangrove_reef_36_vs_30_p0.05.txt",
  "mangrove_reef_36_vs_33_p0.05.txt",
  "reef_reef_33_vs_30_p0.05.txt",
  "reef_reef_36_vs_30_p0.05.txt",
  "reef_reef_36_vs_33_p0.05.txt",
  "wild_mangrove_33_vs_30_p0.05.txt",
  "wild_mangrove_36_vs_30_p0.05.txt",
  "wild_mangrove_36_vs_33_p0.05.txt",
  "wild_reef_33_vs_30_p0.05.txt",
  "wild_reef_36_vs_30_p0.05.txt",
  "wild_reef_36_vs_33_p0.05.txt"
)


for (direction in c("up", "down")) {
  for (file in files) {
    DE <- read.delim(paste("2-Differential_gene_expression/out/pairwise", file, sep = "/"), sep = "\t", header = TRUE)

    if (direction == "up") {
      DE <- filter(DE, stat > 0)
    }
    if (direction == "down") {
      DE <- filter(DE, stat < 0)
    }

    DE$gene_names <- row.names(DE)
    DE_genes <- DE$gene_names
    message("Number of genes in file: ", length(DE_genes))

    # missing <- ! DE_genes %in% geneUniverse
    # missing <- which(missing == TRUE)
    # missing <- DE_genes[missing]

    keep <- DE_genes %in% geneUniverse
    keep <- which(keep == TRUE)
    DE_genes <- DE_genes[keep]
    message("Number of genes in file with GO annotation: ", length(DE_genes))
    # TODO: Look at the ones without annotation and try to find out why

    # make named vector list of factors showing which are GOI
    geneList <- factor(as.integer(geneUniverse %in% DE_genes))
    names(geneList) <- geneUniverse

    # BP
    # create topgo data object
    myGOdata <- new("topGOdata", description = "My project", ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
    # test for significance
    # run weighted algorithm as classic doesnt take into consideration GO hierarchy so could overrepresent enrichment
    resultFisher <- runTest(myGOdata, algorithm = "weight01", statistic = "fisher")
    # generate a table of results using Gentable function
    allGO <- usedGO(object = myGOdata)
    allRes <- GenTable(myGOdata, weightFisher = resultFisher, orderBy = "resultsFisher", ranksOf = "weightFisher", topNodes = length(allGO))
    # change format to non-scientific digits
    options(scipen = 999)
    # correct for multiple testing e.g. a p-adjusted value
    allRes$adjusted.p <- p.adjust(allRes$weightFisher, method = "bonferroni", n = length(allRes$weightFisher))
    allRes$q.value <- p.adjust(allRes$weightFisher, method = "fdr", n = length(allRes$weightFisher))
    BP <- allRes
    BP$ontology <- "BP"

    # MF
    myGOdata <- new("topGOdata", description = "My project", ontology = "MF", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
    resultFisher <- runTest(myGOdata, algorithm = "weight01", statistic = "fisher")
    allGO <- usedGO(object = myGOdata)
    allRes <- GenTable(myGOdata, weightFisher = resultFisher, orderBy = "resultFisher", ranksOf = "weightFisher", topNodes = length(allGO))
    options(scipen = 999)
    allRes$adjusted.p <- p.adjust(allRes$weightFisher, method = "bonferroni", n = length(allRes$weightFisher))
    allRes$q.value <- p.adjust(allRes$weightFisher, method = "fdr", n = length(allRes$weightFisher))
    MF <- allRes
    MF$ontology <- "MF"

    # CC
    myGOdata <- new("topGOdata", description = "My project", ontology = "CC", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
    resultFisher <- runTest(myGOdata, algorithm = "weight01", statistic = "fisher")
    allGO <- usedGO(object = myGOdata)
    allRes <- GenTable(myGOdata, weightFisher = resultFisher, orderBy = "resultFisher", ranksOf = "weightFisher", topNodes = length(allGO))
    options(scipen = 999)
    allRes$adjusted.p <- p.adjust(allRes$weightFisher, method = "bonferroni", n = length(allRes$weightFisher))
    allRes$q.value <- p.adjust(allRes$weightFisher, method = "fdr", n = length(allRes$weightFisher))
    CC <- allRes
    CC$ontology <- "CC"

    out <- rbind(BP, CC, MF)
    out.s <- subset(out, out$weightFisher < 0.001)
    write.table(out.s, paste0("4-Functional_annotation/out/", direction, "_", file), row.names = FALSE, sep = "\t", quote = FALSE)
  }
}


for (direction in c("up", "down")) {
  combined <- NULL
  for (file in files) {
    file_data <- read_tsv(paste0("4-Functional_annotation/out/", direction, "_", file))
    file_data$file <- file
    file_data$regulation <- direction

    combined <- rbind(combined, file_data)
  }
  write_tsv(combined, paste0("4-Functional_annotation/out/combined_", direction, ".tsv"))
}
