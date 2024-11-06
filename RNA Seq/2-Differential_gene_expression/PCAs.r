library(DESeq2)
library(ggplot2)
source("../plotTheme.r")

setwd("~/Documents/PhD/Chapter 4/FINAL DATA AND CODE/RNA Seq")

source("2-Differential_gene_expression/dds_data.r")

levels(as.factor(metadata$Group)) # checks levels in group

design(dds) <- ~TempGroup # identifies comparison of interest
dds <- DESeq(dds)

vst <- varianceStabilizingTransformation(dds)

# preparing the data
pca_data <- plotPCA(vst, intgroup = "TempGroup", returnData = TRUE)
pca_data$SampleID <- metadata$SampleID
pca_data$Temperature <- as.factor(metadata$Temperature)
pca_data$Group <- as.factor(metadata$Group)

percentVar <- round(100 * attr(pca_data, "percentVar"))

# creates plot
ggplot(pca_data, (aes(x = PC1, y = -PC2))) +
  stat_ellipse(geom = "polygon", level=0.8, aes(col = Group, fill = after_scale(alpha(colour, 0.2))), lwd=0) +
  geom_point(aes(shape = Temperature, col = Group, group = Group), size = 3) +
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
  theme(aspect.ratio = 1) + 
  theme(legend.key=element_blank()) +
  plotTheme

# save plot
ggsave("2-Differential_gene_expression/plots/PCA_TempGroupEllipse.png", dpi = 300, width = 9, height = 6)
