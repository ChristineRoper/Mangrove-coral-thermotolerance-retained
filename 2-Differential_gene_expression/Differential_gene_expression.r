setwd("~/Documents/PhD/Chapter 4/Data/R analysis/RNA-seq/P. acuta new genome")

library(DESeq2)
library(ggplot2)
source("../christineTheme.r")

# unload and reload dplyr.
# Because the global environment is shared in R, and libraries overwrite things, something keeps replacing dplyr functions with
# versions that do differnt things, so the code breaks. We need to force them back to the functions we expect.
unloadNamespace("dplyr")
library(dplyr)


source("2-Differential_gene_expression/dds_data.r")

dds <- DESeq(dds)


############### LRT things ###################
output <- NULL

# By Temperature

threshold <- 0
alpha <- 0.05

design(dds) <- ~Temperature
dds <- DESeq(dds, test = "LRT", reduced = ~1)
res <- results(dds)
output <- rbind(output, data.frame(
  by = "Temperature",
  notallzero = sum(res$baseMean > 0),
  up = sum(res$padj < alpha & res$log2FoldChange > threshold, na.rm = TRUE),
  down = sum(res$padj < alpha & res$log2FoldChange < threshold, na.rm = TRUE),
  filt = sum(!is.na(res$pvalue) & is.na(res$padj)),
  outlier = sum(res$baseMean > 0 & is.na(res$pvalue))
))

res %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  select(gene, padj) %>%
  write_csv("2-Differential_gene_expression/out/DEG_by_Temperature.csv")

# By Group

design(dds) <- ~Group
dds <- DESeq(dds, test = "LRT", reduced = ~1)
res <- results(dds)
output <- rbind(output, data.frame(
  by = "Group",
  notallzero = sum(res$baseMean > 0),
  up = sum(res$padj < alpha & res$log2FoldChange > threshold, na.rm = TRUE),
  down = sum(res$padj < alpha & res$log2FoldChange < threshold, na.rm = TRUE),
  filt = sum(!is.na(res$pvalue) & is.na(res$padj)),
  outlier = sum(res$baseMean > 0 & is.na(res$pvalue))
))

res %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  select(gene, padj) %>%
  write_csv("2-Differential_gene_expression/out/DEG_by_Group.csv")

# By Interaction

design(dds) <- ~ Temperature + Group + Temperature:Group
dds <- DESeq(dds, test = "LRT", reduced = ~ Temperature + Group)
res <- results(dds)
output <- rbind(output, data.frame(
  by = "Interaction",
  notallzero = sum(res$baseMean > 0),
  up = sum(res$padj < alpha & res$log2FoldChange > threshold, na.rm = TRUE),
  down = sum(res$padj < alpha & res$log2FoldChange < threshold, na.rm = TRUE),
  filt = sum(!is.na(res$pvalue) & is.na(res$padj)),
  outlier = sum(res$baseMean > 0 & is.na(res$pvalue))
))

res %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  select(gene, padj) %>%
  write_csv("2-Differential_gene_expression/out/DEG_by_Interaction.csv")


write.csv(output, "2-Differential_gene_expression/out/Summary.csv")


############### Differential pairwise comparisons ###################

design(dds) <- ~TempGroup
dds <- DESeq(dds)
resultsNames(dds)

groups <- levels(metadata$Group)
temperatures <- levels(metadata$Temperature)

output <- NULL
for (group in groups) {
  cat("\033[1m", group, "\033[0m\n\n", sep = "")
  for (i in 1:(length(temperatures) - 1)) {
    for (j in (i + 1):length(temperatures)) {
      t2 <- temperatures[i]
      t1 <- temperatures[j]
      # t1 <- temperatures[i]
      # t2 <- temperatures[j]
      cat(t1, "vs", t2, "\n")
      tempGroup1 <- paste(t1, group, sep = "_")
      tempGroup2 <- paste(t2, group, sep = "_")
      # print(paste(tempGroup1, tempGroup2, sep = " vs "))
      res <- results(dds, contrast = c("TempGroup", tempGroup1, tempGroup2))

      threshold <- 0
      alpha <- 0.05
      output <- rbind(output, data.frame(
        a = tempGroup1,
        b = tempGroup2,
        notallzero = sum(res$baseMean > 0),
        up = sum(res$padj < alpha & res$log2FoldChange > threshold, na.rm = TRUE),
        down = sum(res$padj < alpha & res$log2FoldChange < threshold, na.rm = TRUE),
        filt = sum(!is.na(res$pvalue) & is.na(res$padj)),
        outlier = sum(res$baseMean > 0 & is.na(res$pvalue))
      ))

      res %>%
        as.data.frame() %>%
        filter(padj < 0.05) %>%
        write.table(
          paste0("2-Differential_gene_expression/out/pairwise/", group, "_", t1, "_vs_", t2, "_p0.05.txt"),
          sep = "\t",
          quote = FALSE,
          row.names = TRUE
        )

      # res %>%
      #   as.data.frame() %>%
      #   rownames_to_column("gene") %>%
      #   select(gene, padj) %>%
      #   write_csv(paste("GO_MWU/input", group, t1, "vs", t2, ".csv", sep = "_"))

      summary(res, alpha = 0.05)
    }
  }
}

for (temperature in temperatures) {
  cat("\033[1m", temperature, "\033[0m\n\n", sep = "")
  for (i in 1:(length(groups) - 1)) {
    for (j in (i + 1):length(groups)) {
      g1 <- groups[i]
      g2 <- groups[j]
      cat(g1, "vs", g2, "\n")
      tempGroup1 <- paste(temperature, g1, sep = "_")
      tempGroup2 <- paste(temperature, g2, sep = "_")

      res <- results(dds, contrast = c("TempGroup", tempGroup1, tempGroup2))

      threshold <- 0
      alpha <- 0.05
      output <- rbind(output, data.frame(
        a = tempGroup1,
        b = tempGroup2,
        notallzero = sum(res$baseMean > 0),
        up = sum(res$padj < alpha & res$log2FoldChange > threshold, na.rm = TRUE),
        down = sum(res$padj < alpha & res$log2FoldChange < threshold, na.rm = TRUE),
        filt = sum(!is.na(res$pvalue) & is.na(res$padj)),
        outlier = sum(res$baseMean > 0 & is.na(res$pvalue))
      ))

      res %>%
        as.data.frame() %>%
        filter(padj < 0.05) %>%
        write.table(
          paste0("2-Differential_gene_expression/out/pairwise/", temperature, "_", g1, "_vs_", g2, "_p0.05.txt"),
          sep = "\t",
          quote = FALSE,
          row.names = TRUE
        )

      # res %>%
      #   as.data.frame() %>%
      #   rownames_to_column("gene") %>%
      #   select(gene, padj) %>%
      #   write_csv(paste("GO_MWU/input", temperature, g1, "vs", g2, ".csv", sep = "_"))

      summary(res, alpha = 0.05)
    }
  }
}

write.csv(output, "2-Differential_gene_expression/out/Summary_pairwise.csv")

# Bubble plots showing results across pairwise comparisons
output %>% ggplot(aes(x = a, y = b, size = up)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("2-Differential_gene_expression/plots/up2.png")
output %>% ggplot(aes(x = a, y = b, size = down)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("2-Differential_gene_expression/plots/down2.png")


# Scatterplot not very useful
# plotCounts(dds, gene = which.min(res1$padj), intgroup = "Temperature", xlab = "Temperature")
# ggplot(d, aes(x = Temperature, y = count)) +
#   geom_point(position = position_jitter(w = 0.1, h = 0)) +
#   scale_y_log10(breaks = c(200, 400, 800, 1600))
# ggsave('p.png')


####### Up/Down Bars ######

plotData <- read.csv("2-Differential_gene_expression/out/Summary_pairwise.csv") %>%
  tidyr::unite(comparison, c("a", "b"), remove = FALSE) %>%
  separate_wider_delim(a, delim = "_", names = c("temp1", "group1a", "group1b")) %>%
  separate_wider_delim(b, delim = "_", names = c("temp2", "group2a", "group2b")) %>%
  mutate(group1 = paste(group1a, group1b, sep = " to ")) %>%
  mutate(group2 = paste(group2a, group2b, sep = " to ")) %>%
  mutate(groups = paste(group1, group2, sep = " vs ")) %>%
  mutate(temps = paste(temp1, temp2, sep = " vs ")) %>%
  mutate(comparisonType = if_else(temp1 == temp2, "Group", "Temperature")) %>%
  arrange(comparisonType)

p <- plotData %>% ggplot(aes(x = comparison, fill = comparisonType)) +
  geom_bar(aes(weight = up)) +
  geom_bar(aes(weight = -down), alpha = 0.6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("2-Differential_gene_expression/plots/updown2.png", p)

p <- plotData %>%
  filter(comparisonType == "Group") %>%
  ggplot(aes(x = temp1, fill = groups)) +
  geom_bar(aes(weight = up), position = "dodge", col = "#000000", lwd = 0.2) +
  geom_bar(aes(weight = -down), position = "dodge", alpha = 0.4, col = "#000000", lwd = 0.2) +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_discrete(labels = c(
    "mangrove to reef vs reef to reef" = "mangrove to reef vs reef to reef",
    "mangrove to reef vs wild to mangrove" = "mangrove to reef vs wild mangrove",
    "mangrove to reef vs wild to reef" = "mangrove to reef vs wild reef ",
    "reef to reef vs wild to mangrove" = "reef to reef vs wild mangrove ",
    "reef to reef vs wild to reef" = "reef to reef vs wild reef ",
    "wild to mangrove vs wild to reef" = "wild mangrove vs wild reef"
  )) +
  labs(x = "Temperature (ºC)", y = "Count", fill = "Group") +
  christineTheme
ggsave("2-Differential_gene_expression/plots/updown_Group2.png", width = 8, height = 4)

p <- plotData %>%
  filter(comparisonType == "Temperature") %>%
  ggplot(aes(x = group1, fill = temps)) +
  geom_bar(aes(weight = up), position = "dodge", col = "#333333", lwd = 0.3) +
  geom_bar(aes(weight = -down), position = "dodge", alpha = 0.4, col = "#000000", lwd = 0.3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  christineTheme
ggsave("2-Differential_gene_expression/plots/updown_Temperature2.png", width = 5, height = 4)

p <- read.csv("2-Differential_gene_expression/out/Summary.csv") %>%
  ggplot(aes(x = by)) +
  geom_bar(aes(weight = up), fill = "#595959", col = "#000000", lwd = 0.2) +
  geom_bar(aes(weight = -down), fill = "#D7D7D7", col = "#000000", lwd = 0.2) +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1)) + +
  scale_x_discrete(limits = c("Temperature", "Group", "Interaction")) +
  labs(x = "", y = "Count", fill = "Group") +
  christineTheme
ggsave("2-Differential_gene_expression/plots/updown_12.png", width = 3, height = 4)
