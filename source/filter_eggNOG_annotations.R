library(tidyverse)
options(width=110)

# Load data
annotations_5.0 = read_tsv(
  "intermediate/eggNOG.tab",
  col_names=c("Ortholog", "Category", "Description"),
  col_types = "_???",
  quote = ""
) %>% distinct()

annotations_4.1 = read_tsv(
  "intermediate/eggNOG_4.1.tab.gz",
  col_names=c("Ortholog", "Description", "Category"),
  col_types="?__??____"
) %>% distinct()

annotations_4.1 = annotations_4.1 %>%
  mutate(
    Category = sapply(
      str_extract_all(Category, "[A-Z]"),
      paste, collapse=""
    )
  )

annotations_5.0 = mutate(annotations_5.0, Version = "5.0")
annotations_4.1 = mutate(annotations_4.1, Version = "4.1")

orthologs = read_tsv(
  "data/uniprot_eggNOG.tab",
  col_names = c("UniProt_entry", "Ortholog")
) %>%
  mutate(Version = ifelse(str_starts(Ortholog, "ENOG41"), "4.1", "5.0")) %>%
  mutate(Ortholog = str_remove(Ortholog, "^ENOG(50)|^ENOG(41)")) %>%
  select(-UniProt_entry) %>%
  distinct()

# Select relevant annotations and write to file
annotations = bind_rows(
  inner_join(orthologs, annotations_5.0),
  inner_join(orthologs, annotations_4.1)
)

# Save annotations
write_tsv(annotations, "data/eggNOG_annotations.tab")
