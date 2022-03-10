library(here)
library(corncob)
source(here("code/helpful_functions.R"))

## Vibrio plot
da_site_wild <- readRDS(here("data/processed/da_site_wild.rds"))
da_site_wild_clean <- cleanDA(da_site_wild, "site")
# subset only to vibrionales in wild guts
ps_vibrio <- ps %>% subset_samples(sampleType=="gut" & wild_or_F2=="wild") %>% subset_taxa(Order=="Vibrionales") %>% prune_taxa(taxa_sums(.)>0,.)
# subset only to vibrio that appear in more than 3 samples and melt
ps_vibrio_overlap <- ps_vibrio %>% filter_taxa2(~ sum(.>0)>3) %>% merge_samples2(group="site") %>%  transform_sample_counts(~ ./sum(.)) %>% psmelt()
# get significantly different vibrios
sig_vibrios <- da_site_wild_clean %>% filter(Order=="Vibrionales" & measure=="mu") %>% use_series("ASV")
# make colors for different vibrios
cols.vibrio <- c(hcl.colors(length(sig_vibrios), palette="Zissou 1"), gray.colors(44-length(sig_vibrios)))
p_vibrio <- ps_vibrio_overlap %>% 
  ggplot(aes(x="", y=Abundance, fill=OTU)) + 
  geom_bar(width=1, stat="identity") + 
  coord_polar("y", start=0) + 
  facet_wrap(~site) + 
  theme_void() + theme(legend.position = "bottom", legend.title = element_blank(), plot.background = element_rect(fill="white", linetype = 0)) + 
  scale_fill_manual(values=cols.vibrio, breaks=c(sig_vibrios, p$data$OTU %>% unique() %>% .[!. %in% sig_vibrios]))

ggsave(p_vibrio, filename=here("figures/plot_vibrio_pie.png"))
