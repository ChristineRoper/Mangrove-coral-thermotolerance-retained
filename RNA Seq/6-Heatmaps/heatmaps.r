library(readr)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(purrr)

# ensure 'select' is the function from the dplyr package
select <- dplyr::select

setwd("RNA Seq")


### Read and join data ####

# Read the weighted fisher results from functional annotation
Weightfisher_results <- bind_rows(
    read_tsv("4-Functional_annotation/output/combined_up.tsv") %>% mutate(direction = "up"),
    read_tsv("4-Functional_annotation/output/combined_down.tsv") %>% mutate(direction = "down")
) %>%
    select(GO.ID, Term, ontology, weightFisher, direction, file) # keep only the columns we're interested in

# read GO slim files output from previous step
data <- bind_rows(
    read_tsv("5-Goslim/output/up.tsv", col_names = FALSE),
    read_tsv("5-Goslim/output/down.tsv", col_names = FALSE)
) %>%
    # Add column names - the files above have no header row
    rename(slim = X1, gene_count = X2, percentage = X3, slim_description = X4, file = X5, direction = X6, ontology = X7, GOTerms = X8) %>%
    # Separate GO terms into individual rows. They are currently one column (comma separated)
    separate_longer_delim(GOTerms, delim = ", ") %>%
    # add on weighted fisher results
    right_join(Weightfisher_results, by = join_by(GOTerms == GO.ID, ontology == ontology, direction == direction, file == file))

# read in the terms of interest we have identified as most important
terms_of_interest <- read_excel("6-Heatmaps/Interaction term of interest with genes.xlsx", sheet = "Shortened terms") %>%
    select(slim, GOTerms, description_short, christine_description) # keep only the columns we're interested in

selected_data <- data %>%
    # Only keep 36 ºC
    filter(startsWith(file, "group-by-temperature/36_")) %>%
    # Filter out any data that are not in "Interaction term of interest with genes.xlsx"
    semi_join(terms_of_interest, by = c("GOTerms", "slim")) %>%
    # Join data from "Interaction term of interest with genes.xlsx"
    left_join(terms_of_interest, by = c("GOTerms", "slim"))


#### Make the heatmap ####

# find the maximum weighted fisher value - used to determine the colour scale
maximum <- max(selected_data$weightFisher, na.rm = TRUE)

selected_data %>%
    # set the order of comparisons in the heatmap
    mutate(file = factor(file, levels = c(
        "group-by-temperature/36_mangrove_reef_vs_wild_mangrove_p0.05.txt",
        "group-by-temperature/36_mangrove_reef_vs_reef_reef_p0.05.txt",
        "group-by-temperature/36_reef_reef_vs_wild_mangrove_p0.05.txt",
        "group-by-temperature/36_reef_reef_vs_wild_reef_p0.05.txt",
        "group-by-temperature/36_mangrove_reef_vs_wild_reef_p0.05.txt",
        "group-by-temperature/36_wild_mangrove_vs_wild_reef_p0.05.txt"
    ))) %>%
    # We want the lowest significance to be in the center of the scale, with increasing positive significance above and negative below.
    # However, the lowest significance is a higher number (greater p-value). Inverting the values makes the lower significance closer to 0 (center)
    # It is important to add a custom scale to the legend so it shows the real (original) values, not the inverted ones (see line 80 below).
    # We invert each weighted fisher value by subtracting it from the maximum value
    mutate(weightFisher_inv = maximum - weightFisher) %>%
    # We multiply the downregulated values by -1 so they become negative
    # Now the upregulated values are positive and the downregulated values are negative
    mutate(weightFisher_inv = weightFisher_inv * ifelse(direction == "down", -1, 1)) %>%
    # Flip 36_reef_reef_vs_wild_mangrove because other mangrove vs reef comparisons have mangrove first. We will also flip the label below
    mutate(weightFisher_inv = weightFisher_inv * ifelse(file == "group-by-temperature/36_reef_reef_vs_wild_mangrove_p0.05.txt", -1, 1)) %>%
    ggplot(aes(x = file, y = description_short, fill = weightFisher_inv)) +
    geom_tile(color = "gray90", size = 0.05) +
    scale_fill_gradient2(
        low = "#0095ff", mid = "#eeeeee", high = "#ff4400",
        # Range from -maximum to 0 to maximum. Multiply by 0.99 to ensure the labels are just within the legend area
        breaks = c(-maximum * 0.99, 0, maximum * 0.99),
        # Make labels match the original weighted fisher values (not inverted)
        labels = c("< 0.1E-10", maximum, "< 0.1E-10")
    ) +
    scale_y_discrete(position = "right") + # put GO terms on the right side of the heatmap
    facet_grid(
        christine_description ~ file,
        scales = "free",
        space = "free",
        switch = "y",
        labeller = labeller(
            christine_description = label_wrap_gen(width = 30, multi_line = TRUE), # wrap labels larger than 30 characters
            file = c( # rename comparisons for readability
                "group-by-temperature/36_mangrove_reef_vs_reef_reef_p0.05.txt" = "mangrove to reef\nvs\nreef to reef",
                "group-by-temperature/36_mangrove_reef_vs_wild_mangrove_p0.05.txt" = "mangrove to reef\nvs\nwild mangrove",
                "group-by-temperature/36_mangrove_reef_vs_wild_reef_p0.05.txt" = "mangrove to reef\nvs\nwild reef",
                "group-by-temperature/36_reef_reef_vs_wild_mangrove_p0.05.txt" = "wild mangrove\nvs\nreef to reef", # flipped, see above
                "group-by-temperature/36_reef_reef_vs_wild_reef_p0.05.txt" = "reef to reef\nvs\nwild reef",
                "group-by-temperature/36_wild_mangrove_vs_wild_reef_p0.05.txt" = "wild mangrove\nvs\nwild reef"
            )
        )
    ) +
    labs(y = "Go Term", fill = "Weighted Fisher") + # label right side y-axis, fill gives the colouring based on weighted fisher values
    theme(
        strip.text.y.left = element_text(angle = 0),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        strip.background = element_rect(fill = "#cccccc", color = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "#999999", fill = NA, size = unit(0.5, "pt")),
        panel.spacing = unit(-0.5, "pt"),
        axis.text.x = element_blank(),
        legend.position = c(1.68, 0.9),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)
    )

ggsave("6-Heatmaps/output/36C_group_comparisons_heatmap.png", width = 11.2, height = 10, limitsize = FALSE)


### Output CSV ###
# this is outputting the raw data that goes into these heatmaps

go_id_to_name <- read_csv("3-Reference_to_GO/GOid_to_GOterm.csv") # this file maps GO ID's to their descriptive GO terms
# note: the heatmap cut off descriptions past a certain length, this provides the full description

selected_data %>%
    # Join go_id_to_name to get full GO name
    left_join(go_id_to_name, by = join_by(GOTerms == id)) %>%
    # Drop "Term" - the truncted Go name is no use now we have the full one
    select(-Term) %>%
    # Write out
    write_csv("6-Heatmaps/output/heatmap_data_36C_group_comparisons.csv")



#### find genes associated with GO terms #####

# Not used for this analysis
# Needs to be reviewed to ensure the genes are matched to the GO terms correctly

# # gene_to_go maps Gene IDs to GO terms
# gene_to_go <- read_tsv("3-Reference_to_GO/output/gene_to_go.tsv", na = c("", "NA", "-"))
# # columns include gene_id, transcript_id, GOs, Description, Preferred_name

# # Load DEGs
# files <- levels(as.factor(selected_data$file))
# files <- c(files[1])
# degs <- map_df(files, ~ read_tsv(paste0("2-Differential_gene_expression/output/pairwise/", .x)) %>% # Load files from differential gene expression step
#     mutate(file = .x)) %>% # change ".x" column name to "file"
#     left_join(gene_to_go, by = join_by(baseMean == gene_id)) %>% # DEGs do not have GO terms associated with them, so we need to add this for each gene
#     drop_na(GOs) %>% # drop any genes that do not have GO terms
#     select(file, GOs, Description, Preferred_name) %>% # keep just these 4 columns
#     separate_longer_delim(GOs, delim = ",") %>% # Put each gene into it's own row
#     group_by(GOs, file) %>% # If the GO term is ducplicated,
#     summarize(
#         Descriptions = paste(unique(Description), collapse = ","),
#         Preferred_names = paste(unique(Preferred_name), collapse = ",")
#     )

# selected_data %>%
#     # Join go_id_to_name to get full GO name
#     left_join(go_id_to_name, by = join_by(genes == id)) %>%
#     # Drop "Term" - the truncted Go name is no use now we have the full one
#     select(-Term) %>%
#     # Join DEG stuff
#     left_join(degs, by = join_by(file == file, genes == GOs)) %>%
#     # Only unique terms
#     distinct(genes, Preferred_names, .keep_all = TRUE) %>%
#     # Write out
#     write_tsv("5-Goslim/heatmaps/out/selected_data_unique_with_degs_all.tsv")
