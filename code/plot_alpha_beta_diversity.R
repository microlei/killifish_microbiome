library(tidyverse)
library(phyloseq)
library(ggplot2)
library(vegan)
library(cowplot)
library(patchwork)
library(here)

ps <- readRDS(here("data/processed/ps.rds"))
ps.rel <- readRDS(here("data/processed/psrel.rds"))

theme_set(theme_light() + theme(strip.text = element_blank(),
                                axis.text.x.bottom=element_blank(),
                                axis.title.x.bottom = element_blank(),
                                axis.title.y.left = element_text(size=9),
                                legend.title = element_blank(),
                                legend.text = element_text(size=8),
                                legend.position = "none", 
                                plot.title=element_text(size=9, hjust=0.5, face = "bold"),))

p <- ps %>% subset_samples(wild_or_F2=="wild") %>% prune_taxa(taxa_sums(.)>0,.) %>% plot_richness(x = "site", measures=c("Observed"))
w <- ps %>% subset_samples(sampleType=="water") %>% prune_taxa(taxa_sums(.)>0,.) %>%  plot_richness(x = "site", measures=c("Observed"))
rich_plot<- ggplot(p$data, aes(x=site, y=value, color=site)) + geom_boxplot() + geom_point(data=w$data, aes(x=site, y=value), color="black") + 
  facet_wrap(~variable, scales="free_y") + labs(title="Observed taxa", y="") + 
  scale_color_manual(values=c("coral3", "turquoise4"), labels=c("New Bedford Harbor", "Scorton Creek"), name="")

## hacky way to get legend of water
lw <- ggplot(w$data, aes(x=site, y=value, color=site))+geom_point()+scale_color_manual(values=c("black", "black"), labels=c("Water"), breaks=c("New Bedford Harbor"),name="")+theme(legend.position = "top")
legend_w <- get_legend(lw + theme(legend.position = "bottom", legend.margin = margin(), legend.box.margin = margin(), legend.spacing = unit(0,"cm")))

dat <- ps.rel %>% subset_samples(wild_or_F2=="wild" & sampleType=="gut") %>% sample_data() %>% as("data.frame")
disper <- ps.rel %>% subset_samples(wild_or_F2=="wild" & sampleType=="gut") %>% prune_taxa(taxa_sums(.)>0,.) %>% distance("bray") %>% betadisper(d=., dat$site)
df <- data.frame(disper=disper$distances,group=disper$group)
beta_plot <- ggplot(df, aes(x=group, y=disper)) + geom_boxplot(aes(color=group)) + 
  scale_color_manual(values=c("coral3", "turquoise4"), labels=c("New Bedford Harbor", "Scorton Creek"), name="") +
  labs(title="Beta dispersion") + theme(axis.text.y = element_text(), axis.title.y.left=element_blank())
legend <- get_legend(rich_plot + theme(legend.position = "bottom", legend.margin = margin(), legend.box.margin = margin(), legend.spacing = unit(0,"cm")))

patchwork <- (rich_plot + beta_plot)/legend + legend_w + plot_layout(heights=c(1,.05,.05))
patchwork
ggsave(filename = here("figures/plot_diversity.png"), plot = patchwork, width=3, height=3, unit="in")
