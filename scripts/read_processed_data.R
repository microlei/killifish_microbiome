# read all processed data
library(here)

ps <- readRDS(here("data/processed/ps.rds"))
ps.rel <- readRDS(here("data/processed/psrel.rds"))
ps.melted <- readRDS(here("data/processed/psmelted.rds"))
ps.rel.melted <- readRDS(here("data/processed/psrelmelted.rds"))

taxonomy <- read.delim(here("data/processed/taxonomy.txt"), "\t", header=TRUE, row.names=1, check.names = FALSE)
asvs <- read.delim(here("data/processed/ASVs.txt"), "\t", header=TRUE, row.names=1, check.names = FALSE)
metadata <- read.delim(here("data/processed/metadata.csv"), ",", header = TRUE, row.names=1, check.names=FALSE)
