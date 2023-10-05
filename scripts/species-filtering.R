library(tidyverse)

repeats_hits <- read.csv("metadata/binarized_taxid_exceeds_healthy_human_max_20230821.csv")

repeat_hits_filtered <- repeats_hits %>%
  select(-COMP) %>%
  pivot_longer(!Common.Name, names_to = "gene", values_to = "count") %>%
  filter(count > 0)

species_list <- repeat_hits_filtered %>%
  group_by(Common.Name) %>%
  count() %>%
  arrange(desc(n))

species_to_run <- repeat_hits_filtered %>%
  mutate(modf_name = gsub(" ", "_", Common.Name)) %>%
  select(modf_name) %>%
  unique()

write.csv(species_to_run, "metadata/species-to-validate.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
