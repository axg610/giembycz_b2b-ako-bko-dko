library(tidyverse)

runs <- read_tsv("B2B_KO_checksums.txt", col_names = F) %>%
  select(filename=1) %>%
  separate(filename, into = c("samp", "idk", "lane", "run", "idk2"), sep = "_", remove = F) %>%
  na.omit() %>%
  mutate(samp = factor(samp, levels = unique(samp))) %>%
  arrange(samp, lane, run) %>%
  mutate(lanerun = paste(lane, run, sep = "_")) %>%
  select(samp, lanerun, filename) %>%
  pivot_wider(names_from = lanerun, values_from = filename) %>%
  select(-samp) %>%
  write_tsv("pairedEnds.txt")

