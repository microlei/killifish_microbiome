library(tidyverse)
library(phyloseq)
library(SpiecEasi)

ps <- readRDS("../data/processed/ps.rds")

ps.gut.abund <- ps %>% subset_samples(sampleType=="gut") %>% prune_taxa(taxa_sums(.)>50,.)

hpc_params <- list(rep.num = 50, seed = 99999, ncores = 10)

res <- spiec.easi(ps.gut.abund,
                  method = 'glasso',
                  lambda.min.ratio=1e-2,
                  nlambda=20,
                  pulsar.params=hpc_params)

saveRDS(res, file="../data/processed/spieceasi_all.rds")
