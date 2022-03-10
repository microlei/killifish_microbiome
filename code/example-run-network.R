library(tidyverse);library(phyloseq);library(SpiecEasi)

# Load data
load("asv-tables-processed-18102021.RData", verbose = T)

# Load metadata
metadata <- read.delim("samplelist-metadata.txt")

# Select eukaryotes only and create wide format dataframe
insitu_wide_nosingle <- asv_insitu_qc %>%
  filter(Domain == "Eukaryota") %>%
  filter(!grepl("_Plume001_", SAMPLE)) %>% #removing "near vent background", not relevant in other data sets
  select(FeatureID, Taxon, SAMPLE, value) %>%
  pivot_wider(names_from = SAMPLE, values_from = value, values_fill = 0) %>%
  mutate(PREVALENCE = rowSums(select_if(., is.numeric) > 0),
         SEQ_TOTAL = rowSums(select_if(., is.numeric))) %>%
  filter(PREVALENCE >= 1) %>% 
  filter(SEQ_TOTAL >= 100)

# head(insitu_wide)
insitu_samples <- as.character(colnames(insitu_wide_nosingle %>% select(-Taxon, -FeatureID)))


# make matrices for phyloseq
insitu_tax_matrix <- insitu_wide_nosingle %>%
  select(FeatureID, Taxon) %>%
  separate(Taxon, c("Domain", "Supergroup",
                  "Phylum", "Class", "Order",
                  "Family", "Genus", "Species"), sep = ";") %>%
  column_to_rownames(var = "FeatureID") %>%
  as.matrix

insitu_asv_matrix <- insitu_wide_nosingle %>%
  select(-Taxon) %>%
  column_to_rownames(var = "FeatureID") %>%
  as.matrix

# Align row names for each matrix
rownames(insitu_tax_matrix) <- row.names(insitu_asv_matrix)

## Extract relevant metadata information
# head(metadata)
metadata_insitu <- metadata %>%
  filter(SAMPLE %in% insitu_samples) %>% # from reformatting df above
  select(SAMPLE, VENT, SITE, SAMPLETYPE, YEAR) %>%
  unite(SAMPLELABEL, VENT, SITE, SAMPLETYPE, YEAR, sep = "_", remove = FALSE) %>% 
  unite(TYPE_SITE, SITE, SAMPLETYPE, sep = "_", remove = FALSE)

rownames(metadata_insitu) <- metadata_insitu$SAMPLE

# Import asv and tax matrices
ASV = otu_table(insitu_asv_matrix, taxa_are_rows = TRUE)
TAX = tax_table(insitu_tax_matrix)

phylo_obj <- phyloseq(ASV, TAX)

# Import metadata as sample data in phyloseq
samplenames <- sample_data(metadata_insitu)

physeq_insitu = merge_phyloseq(phylo_obj, samplenames)

## Check
physeq_insitu

# Run spiec easi with glasso
pargs2 <- list(rep.num = 50, seed = 10010, ncores = 10)
spec_glasso_microeuk <- spiec.easi(physeq_insitu, method = 'glasso', lambda.min.ratio=1e-2, nlambda=20,pulsar.params=pargs2)

save(spec_glasso_microeuk, file = "spiec-easi-output-03-12-21.RData")

# Extract and add information on weight
#corr_we <- cov2cor(as.matrix(getOptCov(spec_glasso_microeuk))
#colnames(corr_we) <- rownames(corr_we) <- colnames(physeq_insitu)

#weighted_adj_mat <- corr_we*getRefit(spec_glasso_microeuk)
#colnames(weighted_adj_mat) <- rownames(weighted_adj_mat) <- colnames(physeq_insitu)
