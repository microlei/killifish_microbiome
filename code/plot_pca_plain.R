# Code to plot the PCA of all samples and then subplots showing loadings of the top 10 ASVs

library(tidyverse)
library(phyloseq)
library(ggrepel)
library(speedyseq)
library(microbiome)
library(vegan)
library(here)

source(here("code/helpful_functions.R"))

set.seed(100)

# load data
ps <- readRDS(here("data/processed/ps.rds"))
ps.rel <- readRDS(here("data/processed/psrel.rds"))
metadata <- sample_data(ps) %>% as("data.frame")
taxonomy <- ps %>% tax_table %>% as_tibble %>% rename(ASV=".otu")

# Use zComposition to impute zeros before performing a center log ratio transform
f <- zCompositions::cmultRepl(t(otu_table(ps)), method="CZM", label=0, output="p-counts") %>% t()
ps.zcomp <- ps
otu_table(ps.zcomp) <- otu_table(f, taxa_are_rows = TRUE)
ps.zcomp <- transform(ps.zcomp, "clr")

# RDA with no other options is the same as CA/PCA
pca <- rda(t(otu_table(ps.zcomp)))
# Join the sample scores with the metadata
pca.df <- scores(pca, display="sites") %>% data.frame() %>% rownames_to_column(var="sample") %>% left_join(metadata %>% select(sample, site, wild_or_F2, fishType))
# Join the species (ASV) loadings with the taxonomy
pca.species <- scores(pca, display="species") %>% data.frame() %>% rownames_to_column(var="ASV") %>% left_join(taxonomy)
# get summary
pca.summary <- summary(pca)

# normal PCA with no arrows
pca.plot <- ggplot(pca.df, aes(x=PC1, y=PC2)) +
  geom_point(aes(color=fishType), size=4, position="jitter") +
  scale_color_manual(name="", values=cols.fishType) +
  coord_fixed() + theme(legend.position = "none") + 
  labs(x=str_c("PC1 [", round(pca.summary$cont$importance[2,1]*100,2),"%]"),
       y=str_c("PC2 [", round(pca.summary$cont$importance[2,2]*100,2),"%]"))

#ggsave(filename=here("figures/plot_pca_plain.png"), width = 220, height =180, unit="mm" )

# Alternate with ggforce() ellipses
library(ggforce)
pca.plot +
  geom_mark_ellipse(aes(fill=wild_or_F2, label=wild_or_F2), alpha=0, label.fontsize=20) +
  guides(color = guide_legend(ncol=1), fill="none") +
  theme(legend.position = "right", legend.margin = margin(), legend.box.margin = margin(), legend.spacing = unit(0,"cm"), legend.text=element_text(size=11)) +
  coord_cartesian()
ggsave(filename=here("figures/plot_pca_plain.png"))
