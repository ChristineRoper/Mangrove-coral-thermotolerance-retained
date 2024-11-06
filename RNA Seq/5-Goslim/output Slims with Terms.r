library(readr)
library(dplyr)
library(tidyr)

# provides summary dataset for manuscript of GO slims linked with all GO terms enriched by group, temperature, and group at 36 ºC

# ensure 'select' is the function from the dplyr package
select <- dplyr::select

setwd("~/Documents/PhD/Chapter 4/FINAL DATA AND CODE/RNA Seq")


# Read the weighted fisher results from functional annotation
Weightfisher_results <- bind_rows(
    read_tsv("4-Functional_annotation/output/combined_up.tsv") %>% mutate(direction = "up"),
    read_tsv("4-Functional_annotation/output/combined_down.tsv") %>% mutate(direction = "down")
) %>%
    select(GO.ID, Term, ontology, weightFisher, direction, file) # keep only the columns we're interested in

# this file maps GO ID's to their descriptive GO terms
go_id_to_name <- read_csv("3-Reference_to_GO/GOid_to_GOterm.csv")

# read GO slim files output from goslim.r
data <- bind_rows(
    read_tsv("5-Goslim/output/up.tsv", col_names = FALSE),
    read_tsv("5-Goslim/output/down.tsv", col_names = FALSE)
) %>%
    # Add column names - the files above have no header row
    rename(slim = X1, gene_count = X2, percentage = X3, slim_description = X4, file = X5, direction = X6, ontology = X7, GOTerms = X8) %>%
    # Separate GO terms into individual rows. They are currently one column (comma separated)
    separate_longer_delim(GOTerms, delim = ", ") %>%
    # add on weighted fisher results
    right_join(Weightfisher_results, by = join_by(GOTerms == GO.ID, ontology == ontology, direction == direction, file == file)) %>%
    # Join go_id_to_name to get full GO name
    left_join(go_id_to_name, by = join_by(GOTerms == id)) %>%
    # Drop "Term" - the truncted Go name is no use now we have the full one
    select(-Term)


write_csv(
    data %>% filter(startsWith(file, "group-by-temperature/36_")),
    "5-Goslim/output/36C by Group.csv"
)

write_csv(
    data %>% filter(startsWith(file, "by-group/")),
    "5-Goslim/output/by Group.csv"
)

write_csv(
    data %>% filter(startsWith(file, "by-temperature/")),
    "5-Goslim/output/by Temperature.csv"
)
