library(tidyverse)
library(SpiecEasi)
library(igraph)
library(ggraph)
library(tidygraph)
library(microbiome)
library(cowplot)
library(here)
source(here("code/helpful_functions.R"))

ps <- readRDS(here("data/processed/ps.rds"))
ps.wild.abund <- ps %>% subset_samples(sampleType=="gut" & wild_or_F2=="wild") %>% prune_taxa(taxa_sums(.)>50,.)
se_cov <- readRDS("data/processed/se_wild_cov.rds")
se_refit <- readRDS("data/processed/se_wild_refit.rds")

taxonomy <- tax_table(ps.wild.abund) %>% as.data.frame() %>% rownames_to_column(var="ASV")
metadata <- sample_data(ps.wild.abund) %>% as.data.frame()
otus <- otu_table(ps.wild.abund) %>% as.data.frame() %>% rownames_to_column(var="ASV")

colnames(se_cov) <- rownames(se_cov) <- rownames(otu_table(ps.wild.abund))

# get log mean relative abundances of ASVs
ps.rel <- transform_sample_counts(ps, function(x) x/sum(x))
v_abund <- apply(ps.rel %>% subset_taxa(taxa_names(ps.rel) %in% taxa_names(ps.wild.abund)) %>% otu_table(), 1, function(x) log(mean(x)))
# make igraph network
se_net <- adj2igraph(se_refit, vertex.attr = list(name=taxa_names(ps.wild.abund)))

# add vertex attribute
V(se_net)$degree <- degree(se_net)
V(se_net)$abund <- v_abund
V(se_net)$prev <- ps.wild.abund %>% prevalence(detection=0, count=TRUE)
# look at degree distribution
degree.distribution(se_net) %>% data.frame(freq=., deg=0:(length(.)-1)) %>% ggplot(aes(x=deg, y=freq)) + geom_bar(stat="identity")


# Which taxonomic orders act as hubs?
# Which Orders have the highest connectivity (lower right of graphing node degree vs node number)

# Mutate to add node statistics
se_v_flat <- as_tbl_graph(se_net) %>% activate(nodes) %>% mutate(cent_close=centrality_closeness(), cent_hub=centrality_hub(), cent_between = centrality_betweenness(), cent_degree = centrality_degree()) %>% 
  left_join(taxonomy, by=c("name"="ASV"))
order_hub <- se_v_flat %>% as_tibble() %>% group_by(Order) %>% summarise(cent_close=mean(cent_close), cent_hub=mean(cent_hub), cent_between=mean(cent_between))
ggplot(order_hub) + geom_label_repel(aes(x=cent_hub, y=cent_close, size=cent_between, color=cent_between, label=Order))

# Are the most abundant taxa also the ones with the highest hub or centrality score? NO
se_v_flat %>% as_tibble() %>% filter(cent_close>4e-5) %>% pivot_longer(cols=c(cent_close, cent_hub, cent_between, cent_degree), names_to="centrality") %>% 
  ggplot(aes(x=abund, y=value)) + geom_point()+facet_wrap(~centrality, scales="free")

# How does prevalence correlate with centrality scores?
se_v_flat %>% as_tibble() %>% filter(cent_close >4e-5) %>% pivot_longer(cols=c(cent_close, cent_hub, cent_between, cent_degree), names_to="centrality") %>% 
  ggplot(aes(x=prev, y=value)) + geom_point() + facet_wrap(~centrality, scales="free")

# Are there distinct clusters?
se_group <- se_v_flat %>% activate(nodes) %>% mutate(group=group_fast_greedy())
se_group %>% activate(nodes) %>% slice_max(order_by=abund, prop=0.2)
order_hub <- se_v_flat %>% as_tibble() %>% group_by(Order) %>% summarise(nodes = n(), deg = mean(degree), abund = mean(abund))

se_v_flat %>% as_tibble() %>% filter(cent_close>1e-5 & Order=="Vibrionales") %>% ggplot()+geom_label(aes(x=cent_between, y=cent_hub, color=Genus, label=name), position="jitter")

ggplot(se_v_flat)+geom_point(aes(x=abund, y=degree, color=Phylum))

ggplot(order_hub)+geom_text(aes(x=nodes, y=deg, label=Order, size=abund))
# Are the graphs originating from ASV2 different than those originating from other abundant Vibrio ASVs?
# Make ego graphs of ASV2 and other Vibrio ASVs
# Calculate network metrics???



# plot network using igraph
g_layout <- layout_nicely(se_NBH_net)
plot(se_NBH_net, layout=g_layout, vertex.label=NA)



p_net <- plot_network(se_net, ps.wild.abund, type="taxa", color="Phylum", label="value")
p_net_leg <- get_legend(p_net)
p_net +theme(legend.position = "none")
