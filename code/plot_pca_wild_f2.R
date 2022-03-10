# Plots just wild taxa with loadings of the top ten significant differentially abundant taxa between wild and F2 fish
# requres corncob to be run first

library(tidyverse)
library(phyloseq)
library(ggforce)
library(ggrepel)
library(vegan)
library(microbiome)
library(here)

source(here("code/helpful_functions.R"))

set.seed(100)

f <- zCompositions::cmultRepl(t(otu_table(ps)), method="CZM", label=0, output="p-counts") %>% t()
ps.zcomp <- ps
otu_table(ps.zcomp) <- otu_table(f, taxa_are_rows = TRUE)
ps.zcomp <- transform(ps.zcomp, "clr")

metadata <- sample_data(ps.zcomp)%>% as("data.frame") %>% filter(sampleType=="gut")
taxonomy <- tax_table(ps.zcomp) %>% as("data.frame")%>% rownames_to_column(var="ASV")

da_wild_f2 <- readRDS(here("data/processed/da_wild_f2.rds"))
da_wild_f2_clean <- cleanDA(da_wild_f2, "wild_or_F2")

# PCA of just the guts
pca <- rda(t(otu_table(ps.zcomp %>% subset_samples(sampleType=="gut"))))
# Join the sample scores with the metadata
pca.df <- scores(pca, display="sites") %>% data.frame() %>% rownames_to_column(var="sample") %>% left_join(metadata %>% select(sample, site, wild_or_F2, fishType))
# Joing the species (ASV) loadings with the taxonomy
pca.species <- scores(pca, display="species") %>% data.frame() %>% rownames_to_column(var="ASV") %>% left_join(taxonomy)
# get summary
pca.summary <- summary(pca)

# Multiplied PC2 by -1 to get the same orientation as plot_pca_plain.R
pca.plot <- ggplot(pca.df, aes(x=PC1, y=-PC2)) +
  geom_point(aes(color=fishType), size=4, position="jitter") +
  scale_color_manual(name="", values=cols.fishType) +
  coord_fixed() + theme(legend.position = "none") + 
  labs(x=str_c("PC1 [", round(pca.summary$cont$importance[2,1]*100,2),"%]"),
       y=str_c("PC2 [", round(pca.summary$cont$importance[2,2]*100,2),"%]"))

pca.asvs <- da_wild_f2_clean %>% slice_min(order_by = p, n=10)

pca.plot + 
  geom_segment(data=pca.species %>%  filter(ASV %in% pca.asvs$ASV), aes(x=0, y=0, xend=PC1, yend=-PC2), arrow=arrow(angle=25, length=unit(5, "mm")), color="darkgrey") + 
  geom_label_repel(data=pca.species %>% filter(ASV %in% pca.asvs$ASV), aes(x=PC1*1.05, y=-PC2*1.05, label=ASV, fill=Family), alpha=0.75) +
  coord_cartesian() + guides(fill=guide_legend(), color="none") + theme(legend.position = "right", legend.text = element_text(size=20), legend.title = element_text(size=20))

# manually save 900x600 px because ggrepel has a bug
# ggsave(filename=here("figures/pca_wild_f2.png"), width=9, height=6, units="cm")
