library(readr)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(purrr)

setwd("/Users/christineroper/Documents/PhD/Chapter 4/Data/R analysis/RNA-seq/P. acuta new genome")

up <- read_tsv("5-Goslim/out/up.tsv", col_names = FALSE)
down <- read_tsv("5-Goslim/out/down.tsv", col_names = FALSE)

fishers <- bind_rows(
    read_tsv("4-Functional_annotation/out/combined_up.tsv") %>% mutate(direction = "up"),
    read_tsv("4-Functional_annotation/out/combined_down.tsv") %>% mutate(direction = "down")
) %>% select(GO.ID, Term, ontology, weightFisher, direction, file)


data <- bind_rows(up, down) %>%
    rename(slim = X1, gene_count = X2, percentage = X3, slim_description = X4, file = X5, direction = X6, ontology = X7, genes = X8) %>%
    separate_longer_delim(genes, delim = ", ") %>%
    right_join(fishers, by = join_by(genes == GO.ID, ontology == ontology, direction == direction, file == file))


excel_sheets("5-Goslim/heatmaps/Term of interest with genes.xlsx")
sheet <- "Shortened terms"
selected <- read_excel("5-Goslim/heatmaps/Term of interest with genes.xlsx", sheet = sheet) %>%
    dplyr::rename(
        slim = `GO slim`,
        genes = `GO term`,
        description_short = `GO description`
    ) %>%
    select( slim, genes, description_short, christine_description)

cutoff <- 0.0001

selected_data <- data %>%
    filter(startsWith(file, "36_")) %>%
    # Filter out any that are not in "Temps & groups function summary.xlsx"
    semi_join(selected, by = c("genes", "slim")) %>%
    left_join(selected, by = c("genes", "slim"))

# filter(weightFisher < cutoff)

maximum <- max(selected_data$weightFisher, na.rm = TRUE)

max(selected_data$weightFisher)

# levels(as.factor(selected_data$substring))

wrap_long_names_40 <- function(names) {
    str_wrap(names, width = 70)
}

selected_data %>%
    mutate(file = factor(file, levels = c(
        "36_mangrove_reef_vs_wild_mangrove_p0.05.txt",
        "36_mangrove_reef_vs_reef_reef_p0.05.txt",
        "36_reef_reef_vs_wild_mangrove_p0.05.txt",
        "36_reef_reef_vs_wild_reef_p0.05.txt",
        "36_mangrove_reef_vs_wild_reef_p0.05.txt",
        "36_wild_mangrove_vs_wild_reef_p0.05.txt"
    ))) %>%
    # filter(file %in% c(
    #     "36_mangrove_reef_vs_reef_reef_p0.05.txt",
    #     "36_reef_reef_vs_wild_mangrove_p0.05.txt",
    #     "36_mangrove_reef_vs_wild_reef_p0.05.txt",
    #     "36_wild_mangrove_vs_wild_reef_p0.05.txt"
    # )) %>%
    # mutate(slim_description = recode_factor(slim_description,
    #     "carbohydrate metabolic process" = "Carbohydrate metabolic process",
    #     "catalytic activity" = "Catalytic activity",
    #     "chromosome" = "Chromosome",
    #     "DNA binding" = "DNA binding",
    #     "DNA recombination" = "DNA recombination",
    #     "DNA repair" = "DNA repair",
    #     "DNA replication" = "DNA replication",
    #     "generation of precursor metabolites and energy" = "Metabolites and energy",
    #     "immune system process" = "Immune system process",
    #     "lipid metabolic process" = "Lipid metabolic process",
    #     "programmed cell death" = "Programmed cell death",
    #     "protein folding" = "Protein folding",
    #     "signaling" = "Signaling",
    #     "transmembrane transport" = "Transmembrane transport",
    #     "transporter activity" = "Transporter activity",
    #     "wound healing" = "Wound healing"
    # )) %>%
    mutate(weightFisher_inv = maximum - weightFisher) %>%
    # Warning! Inverted now `"down", 1, -1`
    mutate(weightFisher_inv = weightFisher_inv * ifelse(direction == "down", -1, 1)) %>%
    # Flip tih nsone
    mutate(weightFisher_inv = weightFisher_inv * ifelse(file == "36_reef_reef_vs_wild_mangrove_p0.05.txt", -1, 1)) %>%
    # drop_na(name) %>%
    # mutate(gene_direction = paste(gene, direction)) %>%
    ggplot(aes(x = file, y = description_short, fill = weightFisher_inv)) +
    geom_tile(color = "gray90", size = 0.05) +
    scale_fill_gradient2(
        low = "#0095ff", mid = "#eeeeee", high = "#ff4400",
        breaks = c(-maximum * 0.99, 0, maximum * 0.99),
        labels = c("0 (downregulated)", maximum, "0 (upregulated)")
    ) +
    # guides(fill = guide_legend(title = "P-value")) +
    scale_y_discrete(position = "right") + # , labels = wrap_long_names_40) +
    facet_grid(
        christine_description ~ file,
        scales = "free",
        space = "free",
        switch = "y",
        labeller = labeller(
            christine_description = label_wrap_gen(width = 30, multi_line = TRUE),
            # file = c(
            #     "mangrove_reef_30_vs_36_p0.05.txt" = "mangrove-reef 36 vs 30",
            #     "reef_reef_30_vs_36_p0.05.txt" = "reef-reef 36 vs 30",
            #     "wild_mangrove_30_vs_36_p0.05.txt" = "wild mangrove 36 vs 30",
            #     "wild_reef_30_vs_36_p0.05.txt" = "wild reef 36 vs 30"
            # )
            file = c(
                "36_mangrove_reef_vs_reef_reef_p0.05.txt" = "mangrove to reef\nvs\nreef to reef",
                "36_mangrove_reef_vs_wild_mangrove_p0.05.txt" = "mangrove to reef\nvs\nwild mangrove",
                "36_mangrove_reef_vs_wild_reef_p0.05.txt" = "mangrove to reef\nvs\nwild reef",
                "36_reef_reef_vs_wild_mangrove_p0.05.txt" = "wild mangrove\nvs\nreef to reef",
                "36_reef_reef_vs_wild_reef_p0.05.txt" = "reef to reef\nvs\nwild reef",
                "36_wild_mangrove_vs_wild_reef_p0.05.txt" = "wild mangrove\nvs\nwild reef"
            )
        )
    ) +
    labs(y = "Go Term", fill = "Weighted Fisher") +
    theme(
        strip.text.y.left = element_text(angle = 0),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        strip.background = element_rect(fill = "#cccccc", color = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "#999999", fill = NA, size = unit(0.5, "pt")),
        panel.spacing = unit(-0.5, "pt"),
        # strip.clip=False,
        axis.text.x = element_blank(),
        legend.position = c(1.68, 1.065),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)
    )


# ggsave("5-Goslim/heatmaps/Groups_mangrove_vs_reef.png", width = 13, height = 100, limitsize = FALSE)
ggsave("5-Goslim/heatmaps/Groups_mangrove_vs_reef_selected.png", width = 11.2, height = 10, limitsize = FALSE)
    

## ---------

# Output CSV

go_id_to_name <- read_csv("../P. acuta old genome/GOdata/GOname_to_GOid.csv")

selected_data %>%
    # Only unique terms
    distinct(Term, .keep_all = TRUE) %>%
    # Join go_id_to_name to get full GO name
    left_join(go_id_to_name, by = join_by(genes == id)) %>%
    # Drop "Term" - the truncted Go name is no use now we have the full one
    select(-Term) %>%
    # Write out
    write_csv("5-Goslim/heatmaps/out/selected_data_unique_all.csv")


## find DEgenes associated with GO terms

# gene_to_go maps Gene IDs to GO terms
gene_to_go <- read_tsv("3-Reference_to_GO/out/gene_to_go.tsv", na = c("", "NA", "-"))
# gene_id, transcript_id, GOs, Description, Preferred_name

# Load DEGs
files <- levels(as.factor(selected_data$file))
degs <- map_df(files, ~ read_tsv(paste0("2-Differential_gene_expression/out/pairwise/", .x)) %>%
    mutate(file = .x)) %>%
    left_join(gene_to_go, by = join_by(baseMean == gene_id)) %>%
    drop_na(GOs) %>%
    select(file, GOs, Description, Preferred_name) %>%
    separate_longer_delim(GOs, delim = ",") %>%
    group_by(GOs, file) %>%
    summarize(
        Descriptions = paste(unique(Description), collapse = ","),
        Preferred_names = paste(unique(Preferred_name), collapse = ",")
    )

selected_data %>%
    # Join go_id_to_name to get full GO name
    left_join(go_id_to_name, by = join_by(genes == id)) %>%
    # Drop "Term" - the truncted Go name is no use now we have the full one
    select(-Term) %>%
    # Join DEG stuff
    left_join(degs, by = join_by(file == file, genes == GOs)) %>%
    # Only unique terms
    distinct(genes, Preferred_names, .keep_all = TRUE) %>%
    # Write out
    write_tsv("5-Goslim/heatmaps/out/selected_data_unique_with_degs_all.tsv")
