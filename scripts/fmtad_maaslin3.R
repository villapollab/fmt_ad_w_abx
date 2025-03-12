#fmtad_wabx

library(maaslin3)
library(dplyr)
library(phyloseq)

# Load your phyloseq object
ps <- readRDS("data/fmtad_wabx_physeq.rds")

# ----- Create the Feature Table for Maaslin3 -----

# Extract the OTU table as a data frame
otu_df <- as.data.frame(otu_table(ps))

# If taxa are stored as columns, transpose the table so rows correspond to OTUs
if (!taxa_are_rows(ps)) {
  otu_df <- as.data.frame(t(otu_df))
}

# Add the OTU IDs (rownames) as a column
otu_df <- cbind(otu = rownames(otu_df), otu_df)

# Convert the phyloseq taxonomy table to a data frame
tax_df <- as.data.frame(tax_table(ps))

# Check that the "Genus" column exists
if (!"Genus" %in% colnames(tax_df)) {
  stop("The taxonomy table does not contain a column named 'Genus'. Please verify your phyloseq object.")
}

# Add Genus information for each OTU by matching OTU IDs
otu_df$Genus <- tax_df[otu_df$otu, "Genus"]

# Rearrange columns to have "Genus" first (dropping the original OTU column)
otu_df <- otu_df[, c("Genus", setdiff(colnames(otu_df), c("Genus", "otu")))]

# Aggregate (sum) counts by Genus so each Genus appears only once
otu_df_agg <- otu_df %>%
  group_by(Genus) %>%
  summarise(across(everything(), sum)) %>%
  as.data.frame()

# Set rownames to the Genus names and remove the Genus column
rownames(otu_df_agg) <- otu_df_agg$Genus
otu_df_agg$Genus <- NULL

# Write the aggregated feature table with rownames
write.table(otu_df_agg,
            file = "data/maaslin3_feature_table_updated.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE,
            col.names = NA)

# ----- Create the Metadata File -----

# Extract sample metadata
metadata_df <- as.data.frame(sample_data(ps))
metadata_df <- cbind(sample = rownames(metadata_df), metadata_df)

# Write the metadata to a TSV file
write.table(metadata_df,
            file = "data/maaslin3_metadata.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)

# ----- Read Files Back In (Ensuring Rownames are Preserved) -----

taxa_table <- read.delim("data/maaslin3_feature_table_updated.tsv", 
                         check.names = FALSE,
                         header = TRUE,
                         row.names = 1)

metadata <- read.delim("data/maaslin3_metadata.tsv", 
                       check.names = FALSE,
                       header = TRUE,
                       row.names = 1)

metadata$Timepoint <-
  factor(metadata$Timepoint, levels = c('BA', '1DPI', '3DPI', '11DPI'))

# Fit models
fit_out <- maaslin3(input_data = taxa_table,
                    input_metadata = metadata,
                    output = 'data/fmtad_wabx_maaslin3',
                    formula = 'Timepoint',
                    normalization = 'TSS',
                    transform = 'LOG',
                    augment = TRUE,
                    standardize = TRUE,
                    max_significance = 0.1,
                    median_comparison_abundance = TRUE,
                    median_comparison_prevalence = FALSE,
                    max_pngs = 100,
                    cores = 1,
                    save_models = TRUE)


###############################################