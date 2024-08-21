library(DESeq2)
library(ggplot2)

setwd("~/Documents/PhD/Chapter 4/Data/R analysis/RNA-seq/P. acuta new genome")

source("2-Differential_gene_expression/dds_data.r")

levels(as.factor(metadata$Group))

pca <- function(dds, intgroup) {
  vst <- varianceStabilizingTransformation(dds)

  pca_data <- plotPCA(vst, intgroup = intgroup, returnData = TRUE)
  pca_data$SampleID <- metadata$SampleID
  pca_data$Temperature <- as.factor(metadata$Temperature)
  pca_data$Group <- as.factor(metadata$Group)

  percentVar <- round(100 * attr(pca_data, "percentVar"))

  ggplot(pca_data, (aes(x = PC1, y = -PC2))) +
    geom_point(aes(shape = Temperature, col = Group), size = 3) +
    theme_classic() +
    scale_color_discrete(labels = c(
      "mangrove_reef" = "mangrove-reef",
      "reef_reef" = "reef-reef",
      "wild_mangrove" = "wild mangrove",
      "wild_reef" = "wild reef"
    )) +
    scale_shape_discrete(labels = c(
      "30" = "30 ºC",
      "33" = "33 ºC",
      "36" = "36 ºC"
    )) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    theme(axis.text = element_text(size = 20)) +
    theme(axis.title = element_text(size = 20)) +
    theme(legend.title = element_text(size = 20)) +
    theme(text = element_text(size = 20)) +
    theme(aspect.ratio = 1)
}

design(dds) <- ~TempGroup
dds <- DESeq(dds)
pca(dds, "TempGroup") + christineTheme #+ theme(aspect.ratio = 3 / 4)
ggsave("2-Differential_gene_expression/plots/PCA_TempGroup.png", dpi = 300, width = 9, height = 6)

design(dds) <- ~Temperature
dds <- DESeq(dds)
pca(dds, "Temperature")
ggsave("2-Differential_gene_expression/plots/PCA_Temperature.png", dpi = 300, width = 8, height = 8)

design(dds) <- ~Group
dds <- DESeq(dds)
pca(dds, "Group")
ggsave("2-Differential_gene_expression/plots/PCA_Group.png", dpi = 300, width = 8, height = 8)
