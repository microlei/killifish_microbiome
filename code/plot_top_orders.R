# Code to plot the top 5 most abundant orders

library(tidyverse)
library(magrittr)
library(ggplot)
library(here)

source(here("code/helpful_functions.R"))

ps.melted <- readRDS(here("data/processed/psmelted.rds"))
ps <- readRDS(here("data/processed/ps.rds"))

# in psmelt, group by Order and fishType and summarise average relative abundance, taking the top 5 from each fish type
top5_order <- ps.melted %>% filter(sampleType=="gut") %>% group_by(Order, fishType) %>% 
  summarise(mean_relAbund = mean(abund_rel)) %>% group_by(fishType) %>% 
  top_n(5, mean_relAbund) %>% pull(Order)

# glom ps by order and transform to relative abundance, then melting again
ps_top5order_melt <- ps %>% tax_glom("Order") %>% transform_sample_counts(function(x) x/sum(x)) %>% psmelt() %>% 
  filter(Order %in% top5_order & sampleType=="gut")

# plotting all 15 top Orders, observed that some orders have a high average just because of a few outliers
ps_top5order_melt %>% ggplot(aes(x=fishType, y=Abundance)) + geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color=fishType), height = 0, width=.2) + facet_wrap(~Order, scales="free")

# better plot with the clear outlier orders removed
ps_top5order_melt$fishType %<>% factor(factor.fishType) 
ps_top5order_melt %>% filter(Order %in% c("Vibrionales", "Mycoplasmatales", "Lactobacillales", "Clostridiales", "Cyanobacteriales")) %>% 
  ggplot(aes(x=fishType, y=Abundance)) + geom_boxplot(outlier.shape = NA) + geom_jitter(aes(color=fishType), height = 0, width=.2) + facet_wrap(~Order, scales="free") + 
  scale_color_manual(name="Fish Type", values = cols.fishType[1:4]) + theme(axis.text.x = element_blank()) + labs(y="Relative abundance")

ggsave(filename="figures/plot_top_orders.png")
