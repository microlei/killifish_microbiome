ASV differences between wild and F2
================
Lei Ma
Last compiled on 01 October, 2021

``` r
ps <- readRDS(here("data/processed/ps.rds"))
ps.rel <- readRDS(here("data/processed/psrel.rds"))
# ps.rel.melted <- readRDS(here("data/processed/psrelmelted.rds"))
```

## Overlap of ASVs

I was curious how many ASVs were unique to each environment/fish type:
Scorton Creek wild, New Bedford Harbor wild, and the respective F2 fish

``` r
taxa_list <- list(sc_w = subset_samples(ps, fishType=="Scorton Creek wild") %>% prune_taxa(taxa_sums(.)>0,.) %>% taxa_names(),
                  nb_w = subset_samples(ps, fishType=="New Bedford Harbor wild") %>% prune_taxa(taxa_sums(.)>0,.) %>% taxa_names(),
                  sc_f2 = subset_samples(ps, fishType=="Scorton Creek F2") %>% prune_taxa(taxa_sums(.)>0,.) %>% taxa_names(),
                  nb_f2 = subset_samples(ps, fishType=="New Bedford Harbor F2") %>% prune_taxa(taxa_sums(.)>0,.) %>% taxa_names())

taxa_overlap <- calculate.overlap(taxa_list)

venn.diagram(taxa_list, filename = here("figures/Venn_gut.png"), imagetype="png", main="Venn diagram of taxa shared between fish types", category.names = c("SC wild", "NBH wild", "SC F2", "NBH F2"))
```

    ## [1] 1

``` r
knitr::include_graphics(path=here("figures/Venn_gut.png"))
```

<img src="/Users/Lei/Dropbox (MIT)/Apprill Lab/Killifish/killifish_microbiome/figures/Venn_gut.png" width="3000" />

## Core taxa - present in all fish types

Core taxa in wild type are either Vibrionaceae or Mycoplasma while core
taxa in F2 are Enterococcus

``` r
ps.core <- prune_taxa(taxa_names(ps.rel) %in% taxa_overlap$a6, ps.rel) %>% subset_samples(sampleType=="gut")
core.plot <- ps.core %>% plot_bar(fill="Order")+ geom_bar(aes(fill=Order), stat="identity", position="stack")+facet_wrap(~fishType, scales="free_x")
core.plot.legend <- get_legend(core.plot+ guides(fill = guide_legend(nrow = 3))+theme(legend.position = "bottom"))
p <- core.plot + theme(legend.position = "none")
plot_grid(p, core.plot.legend, ncol=1, rel_heights = c(1,.2))
```

![](3_ASV_differences_files/figure-gfm/core-1.png)<!-- -->

``` r
ps.subset <- subset_taxa(ps.rel, Order %in% c("Vibrionales", "Lactobacillales","Mycoplasmatales")) %>% subset_samples(sampleType=="gut")
ps.subset %>% plot_bar(fill="Order")+ geom_bar(aes(fill=Order), stat="identity", position="stack")+facet_wrap(~fishType, scales="free_x")
```

![](3_ASV_differences_files/figure-gfm/vibrio%20myco%20and%20entero-1.png)<!-- -->

Are the vibrios/mycoplasmatales and enterococcus negative associated
with each other within wild types and F2s? Or is it just that
enterococcus overtook the vibrios? Analysis TBD.
