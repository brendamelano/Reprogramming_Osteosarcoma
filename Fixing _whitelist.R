library(dplyr)

whitelist <- read.csv("Perturb_seq_analysis/Perturb_osteo_whitelist_guides.csv")


whitelist <- whitelist %>%
  group_by(target_gene_name) %>%
  mutate(
    id = case_when(
      row_number() == 1 ~ paste0(id, '-1i'),
      row_number() == 2 ~ paste0(id, '-2i'),
      TRUE ~ as.character(id)
    ),
    name = case_when(
      row_number() == 1 ~ paste0(name, '-1i'),
      row_number() == 2 ~ paste0(name, '-2i'),
      TRUE ~ as.character(name)
    )
  ) %>%
  ungroup()

write.csv(whitelist, '~/Desktop/Reprogramming_Osteosarcoma/Perturb_seq_analysis/Perturb_osteo_whitelist_guides.csv')

