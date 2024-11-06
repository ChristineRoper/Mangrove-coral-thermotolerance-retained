library(readr)
library(stringr)

transcript_lengths <- read_tsv(
    "reference genome/Pocillopora_acuta_HIv2.genes.gtf",
    skip = 1, # skip row one, contains no data
    col_names = c( # columns do not have names so here we are giving them names
        "SequenceID", "Source",
        "FeatureType", "FeatureStart", "FeatureEnd",
        "Score", "Stand", "Phase",
        "Attributes"
    )
) %>%
    # Only keep the "transcripts"
    filter(FeatureType == "transcript") %>%
    # Split it to get the gene_id and transcript_id, if there are more attributes, drop them.
    separate_wider_delim(Attributes, delim = regex("; | "), names = c(NA, "gene_id", NA, "transcript_id"), too_many = "drop") %>%
    # Strip the quotemarks "Pocillopora_acuta_HIv2___TS.g24026.t1" => Pocillopora_acuta_HIv2___TS.g24026.t1
    mutate(gene_id = gsub('"', "", gene_id)) %>%
    # Calculate length. The +1 is added because the length should include both endpoints.
    mutate(length = FeatureEnd - FeatureStart + 1) %>%
    # we only need the ID and length
    select(gene_id, length) %>%
    column_to_rownames(var = "gene_id")

raw_counts <- read_csv("1-Process-HTSeq/output/HTSeqCounts.csv") %>% column_to_rownames(var = "Gene")

# normalise raw counts to transcript_lengths
normalised_counts <- data.matrix(raw_counts) / transcript_lengths$length

# Convert to TPM
tpm <- t(t(normalised_counts) * 1e6 / colSums(normalised_counts))

tpm %>%
    as.data.frame() %>%
    rownames_to_column('transcript_id')%>%
    write_csv("1-Process-HTSeq/output/TPM.csv")
