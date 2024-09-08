setwd("~/Documents/PhD/Chapter 4/FINAL DATA AND CODE/RNA Seq")

library(DESeq2)
library(ggplot2)
source("../christineTheme.r")

# unload and reload dplyr.
# Because the global environment is shared in R, and libraries overwrite things, something keeps replacing dplyr functions with
# versions that do differnt things, so the code breaks. We need to force them back to the functions we expect.
unloadNamespace("dplyr")
library(dplyr)

# run dds_data.r which we used to create a DESeq dataset from our data
source("2-Differential_gene_expression/dds_data.r")

# dds <- DESeq(dds)


############### Likelihood ratio test ###################
## Use this to compare factors with 3 or more levels ##
# Note: the LRT tests were only used to summarise up/downregulated genes by temp, group and their interaction, summarised in a figure
output <- NULL

# By Temperature

threshold <- 0 # determine up/downregulated changes - see line 36
alpha <- 0.05 # the significance level for the log2FoldChange

design(dds) <- ~Temperature # set design to be by temperature
dds <- DESeq(dds, test = "LRT", reduced = ~1) # run likelihood ratio test
res <- results(dds)
# put results into output
output <- rbind(output, data.frame(
  by = "Temperature",
  notallzero = sum(res$baseMean > 0),
  up = sum(res$padj < alpha & res$log2FoldChange > threshold, na.rm = TRUE), # if the log2FoldChange is > 0 it's upregulated, if its < 0 it's downregulated
  down = sum(res$padj < alpha & res$log2FoldChange < threshold, na.rm = TRUE),
  filt = sum(!is.na(res$pvalue) & is.na(res$padj)),
  outlier = sum(res$baseMean > 0 & is.na(res$pvalue))
))

res %>%
  as.data.frame() %>% # turn results into dataframe
  rownames_to_column("gene") %>% # turn row names into column called "gene"
  select(gene, padj) %>%
  write_csv("2-Differential_gene_expression/output/DEG_by_Temperature.csv") # write to csv

# By Group

design(dds) <- ~Group
dds <- DESeq(dds, test = "LRT", reduced = ~1)
res <- results(dds)
# add the group results to the temperature results from above
output <- rbind(output, data.frame(
  by = "Group",
  notallzero = sum(res$baseMean > 0),
  up = sum(res$padj < alpha & res$log2FoldChange > threshold, na.rm = TRUE),
  down = sum(res$padj < alpha & res$log2FoldChange < threshold, na.rm = TRUE),
  filt = sum(!is.na(res$pvalue) & is.na(res$padj)),
  outlier = sum(res$baseMean > 0 & is.na(res$pvalue))
))

res %>% # as above
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  select(gene, padj) %>%
  write_csv("2-Differential_gene_expression/output/DEG_by_Group.csv")

# By Interaction

design(dds) <- ~ Temperature + Group + Temperature:Group
dds <- DESeq(dds, test = "LRT", reduced = ~ Temperature + Group)
res <- results(dds)
# add interaction results to above results
output <- rbind(output, data.frame(
  by = "Interaction",
  notallzero = sum(res$baseMean > 0),
  up = sum(res$padj < alpha & res$log2FoldChange > threshold, na.rm = TRUE),
  down = sum(res$padj < alpha & res$log2FoldChange < threshold, na.rm = TRUE),
  filt = sum(!is.na(res$pvalue) & is.na(res$padj)),
  outlier = sum(res$baseMean > 0 & is.na(res$pvalue))
))

res %>% # As above
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  select(gene, padj) %>%
  write_csv("2-Differential_gene_expression/output/DEG_by_Interaction.csv")

# create csv which summarises all LTR results (temperature, group and their interaction) which is later visualised in a figure
write.csv(output, "2-Differential_gene_expression/output/Summary_LRT.csv") # output is the summary of gene counts for each factor


############### Differential pairwise comparisons ###################

groups <- levels(metadata$Group)
temperatures <- levels(metadata$Temperature)

output <- NULL

# Compare by group

design(dds) <- ~Group
dds <- DESeq(dds)
resultsNames(dds) # prints the names its using for each group


# this forloop is making "i" and "j" for group comparisons
for (i in 1:(length(groups) - 1)) {
  for (j in (i + 1):length(groups)) {
    g1 <- groups[i] # group 1
    g2 <- groups[j] # group 2
    message(g1, " vs ", g2)

    res <- results(dds, contrast = c("Group", g1, g2)) # for pairwise comparisons, you can only compare two groups at a time

    # puts results from each pairwise comparison into the dataframe
    threshold <- 0 # as above
    alpha <- 0.05
    output <- rbind(output, data.frame(
      a = g1,
      b = g2,
      notallzero = sum(res$baseMean > 0),
      up = sum(res$padj < alpha & res$log2FoldChange > threshold, na.rm = TRUE),
      down = sum(res$padj < alpha & res$log2FoldChange < threshold, na.rm = TRUE),
      filt = sum(!is.na(res$pvalue) & is.na(res$padj)),
      outlier = sum(res$baseMean > 0 & is.na(res$pvalue))
    ))

    res %>%
      as.data.frame() %>%
      filter(padj < 0.05) %>% # only include padjust < 0.05 in the saved txt file
      write.table(
        paste0("2-Differential_gene_expression/output/pairwise/by-group/", g1, "_vs_", g2, "_p0.05.txt"),
        sep = "\t",
        quote = FALSE,
        row.names = TRUE
      )

    summary(res, alpha = 0.05)
  }
}

# Compare by temperature

design(dds) <- ~Temperature
dds <- DESeq(dds)
resultsNames(dds)

# notes as above
for (i in 1:(length(temperatures) - 1)) {
  for (j in (i + 1):length(temperatures)) {
    t2 <- temperatures[i]
    t1 <- temperatures[j]
    message(t1, " vs ", t2)
    res <- results(dds, contrast = c("Temperature", t1, t2))

    threshold <- 0
    alpha <- 0.05
    output <- rbind(output, data.frame(
      a = t1,
      b = t2,
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
        paste0("2-Differential_gene_expression/output/pairwise/by-temperature/", t1, "_vs_", t2, "_p0.05.txt"),
        sep = "\t",
        quote = FALSE,
        row.names = TRUE
      )

    summary(res, alpha = 0.05)
  }
}

# write.csv(output, "2-Differential_gene_expression/output/Summary_pairwise2.csv")


# Compare between groups within each temperature

design(dds) <- ~TempGroup
dds <- DESeq(dds)
resultsNames(dds)

for (temperature in temperatures) {
  message(temperature)
  for (i in 1:(length(groups) - 1)) {
    for (j in (i + 1):length(groups)) {
      g1 <- groups[i]
      g2 <- groups[j]
      message(g1, " vs ", g2)
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
          paste0("2-Differential_gene_expression/output/pairwise/group-by-temperature/", temperature, "_", g1, "_vs_", g2, "_p0.05.txt"),
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

# # Compare between temperatures within each group - this analysis not used in manuscript

# design(dds) <- ~TempGroup
# dds <- DESeq(dds)
# resultsNames(dds)

# for (group in groups) {
#   message(group)
#   for (i in 1:(length(temperatures) - 1)) {
#     for (j in (i + 1):length(temperatures)) {
#       t2 <- temperatures[i]
#       t1 <- temperatures[j]
#       message(t1, " vs ", t2)
#       tempGroup1 <- paste(t1, group, sep = "_")
#       tempGroup2 <- paste(t2, group, sep = "_")
#       res <- results(dds, contrast = c("TempGroup", tempGroup1, tempGroup2))

#       threshold <- 0 # as above
#       alpha <- 0.05
#       output <- rbind(output, data.frame(
#         a = tempGroup1,
#         b = tempGroup2,
#         notallzero = sum(res$baseMean > 0),
#         up = sum(res$padj < alpha & res$log2FoldChange > threshold, na.rm = TRUE),
#         down = sum(res$padj < alpha & res$log2FoldChange < threshold, na.rm = TRUE),
#         filt = sum(!is.na(res$pvalue) & is.na(res$padj)),
#         outlier = sum(res$baseMean > 0 & is.na(res$pvalue))
#       ))

#       res %>%
#         as.data.frame() %>%
#         filter(padj < 0.05) %>%
#         write.table(
#           paste0("2-Differential_gene_expression/output/pairwise/", group, "_", t1, "_vs_", t2, "_p0.05.txt"),
#           sep = "\t",
#           quote = FALSE,
#           row.names = TRUE
#         )

#       # res %>%
#       #   as.data.frame() %>%
#       #   rownames_to_column("gene") %>%
#       #   select(gene, padj) %>%
#       #   write_csv(paste("GO_MWU/input", group, t1, "vs", t2, ".csv", sep = "_"))

#       summary(res, alpha = 0.05)
#     }
#   }
# }

# outputs results from all pairwise comparisons in csv
write.csv(output, "2-Differential_gene_expression/output/Summary_pairwise.csv")



####### up/down bar plots ######

plotData <- read.csv("2-Differential_gene_expression/output/Summary_pairwise.csv") %>%
  filter(startsWith(a, "30_") | startsWith(a, "33_") | startsWith(a, "36_")) %>%
  tidyr::unite(comparison, c("a", "b"), remove = FALSE) %>%
  separate_wider_delim(a, delim = "_", names = c("temp1", "group1a", "group1b")) %>%
  separate_wider_delim(b, delim = "_", names = c("temp2", "group2a", "group2b")) %>%
  mutate(group1 = paste(group1a, group1b, sep = " to ")) %>%
  mutate(group2 = paste(group2a, group2b, sep = " to ")) %>%
  mutate(groups = paste(group1, group2, sep = " vs ")) %>%
  mutate(temps = paste(temp1, temp2, sep = " vs ")) %>%
  mutate(comparisonType = if_else(temp1 == temp2, "Group", "Temperature")) %>%
  arrange(comparisonType)

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
ggsave("2-Differential_gene_expression/plots/updown_plot_by_temp.png", width = 8, height = 4)


p <- read.csv("2-Differential_gene_expression/output/Summary_LRT.csv") %>% # uses csv saved during LRT (not pairwise analysis)
  ggplot(aes(x = by)) +
  geom_bar(aes(weight = up), fill = "#595959", col = "#000000", lwd = 0.2) +
  geom_bar(aes(weight = -down), fill = "#D7D7D7", col = "#000000", lwd = 0.2) +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1)) + +
  scale_x_discrete(limits = c("Temperature", "Group", "Interaction")) +
  labs(x = "", y = "Count", fill = "Group") +
  christineTheme
ggsave("2-Differential_gene_expression/plots/updown_plot_temp_group_interaction.png", width = 3, height = 4)

