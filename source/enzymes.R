options(width=110)
library(tidyverse)


# Load data
lipsmap = read_tsv(
  "data/annotated_comparison_results.2021-05-17.tab.gz",
  col_types = cols(Protein.names = col_character())
) %>% rename(Locus = Gene_names_ordered_locus)

# Save UniProt IDs of peptides with missing loci
write(
  filter(lipsmap, is.na(Locus)) %>% pull(UniProt_entry) %>% unique(),
  "data/missing_locus_uniprot_IDs.txt"
)

# Load missing locus IDs
uniprot_locus = read_tsv(
  "data/uniprot_locus.tab",
  col_names = c("UniProt_entry", "Locus")
)

# Repair locus ID mapping (for Synechocystis)
lipsmap = bind_rows(
  lipsmap %>%
    filter(UniProt_entry %in% uniprot_locus$UniProt_entry) %>%
    select(-Locus) %>%
    left_join(uniprot_locus),
  lipsmap %>%
    filter(!(UniProt_entry %in% uniprot_locus$UniProt_entry))
)

# Determine the number of detected genes
detected_genes = lipsmap %>%
  group_by(Organism) %>%
  summarise(Detected = length(unique(Locus)))

detected_per_met = lipsmap %>%
  group_by(Organism, Metabolite) %>%
  summarise(Detected = length(unique(Locus)))


# Compare number of significant interactions in EC and non-EC genes
gene_interactions = lipsmap %>%
  select(Organism, Metabolite, Peptide_ID, Locus, EC_number, Sign) %>%
  mutate(Significant = ifelse(Sign == "sign", 1, 0)) %>%
  select(-Sign) %>%
  group_by(Organism, Metabolite, Locus, EC_number) %>%
  summarise(
    Peptides = length(Significant),
    Significant = sum(Significant)
  ) %>%
  mutate(
    Enzyme = factor(
      ifelse(is.na(EC_number), "Other", "Enzyme"),
      levels=c("Enzyme", "Other")
    ),
    Interaction = factor(
      ifelse(Significant > 0, "Interaction", "No interaction"),
      levels=c("Interaction", "No interaction")
    )
  )

# Determine combinations of Organism and Metabolite
org_met = gene_interactions %>%
  ungroup() %>%
  select(Organism, Metabolite) %>%
  distinct()

# Perform Fisher's Exact Test for each metabolite in each organism
fisher_results = bind_rows(lapply(
  1:nrow(org_met),
  function (i) {
    # Select Organism and Metabolite for test
    org_met_test = org_met[i,]
    test_tb = gene_interactions %>%
      filter(
        Organism == org_met_test$Organism,
        Metabolite == org_met_test$Metabolite
      )
    # Calculate contingency table
    c_table = table(test_tb$Enzyme, test_tb$Interaction) %>% as.matrix
    # Calculate interaction fraction
    frac_enzyme = c_table["Enzyme","Interaction"] / sum(c_table["Enzyme",])
    frac_other  = c_table["Other","Interaction"] / sum(c_table["Other",])
    # Perform test
    p = fisher.test(c_table)$p.value
    # Return results
    tibble(
      Organism = org_met_test$Organism,
      Metabolite = org_met_test$Metabolite,
      Enzyme = frac_enzyme,
      Other = frac_other,
      p = p
    )
  }
))

# Adjust p-values
fisher_results = fisher_results %>%
  mutate(q = p.adjust(p, "BH")) %>%
  arrange(q)

# Save results table
write_tsv(
  fisher_results, "results/Fisher_exact_test_for_enzyme_interactions.tab"
)
