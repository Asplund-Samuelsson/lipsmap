options(width=110)
library(tidyverse)

# Load data
lipsmap = read_tsv("data/annotated_comparison_results.tab.gz")

kegg_uniprot = read_tsv("data/KEGGgene_uniprot_organism.tab") %>%
  # Clean up the UniProt ID
  mutate(UniProt_entry = str_remove(UniProt, "up:"))

kegg_ec = read_tsv("data/KEGGgene_EC_organism.tab")

# Determine the number of detected genes
detected_genes = lipsmap %>%
  group_by(Organism) %>%
  summarise(Detected = length(unique(Locus)))

detected_per_met = lipsmap %>%
  group_by(Organism, Metabolite) %>%
  summarise(Detected = length(unique(Locus)))

# Make table of enzymes
enzymes = kegg_uniprot %>% inner_join(kegg_ec)

# Compare number of significant interactions in EC and non-EC genes
gene_interactions = lipsmap %>%
  select(Organism, Metabolite, Peptide_ID, UniProt_entry, Conc, Sign) %>%
  mutate(Significant = ifelse(Sign == "sign", 1, 0)) %>%
  select(-Sign) %>%
  group_by(Organism, Metabolite, Conc, UniProt_entry) %>%
  summarise(
    Peptides = length(Significant),
    Significant = sum(Significant)
  ) %>%
  mutate(
    Enzyme = factor(
      ifelse(UniProt_entry %in% enzymes$UniProt_entry, "Enzyme", "Other"),
      levels=c("Enzyme", "Other")
    ),
    Interaction = factor(
      ifelse(Significant > 0, "Interaction", "No interaction"),
      levels=c("Interaction", "No interaction")
    )
  )

# Determine combinations of Organism, Metabolite, and Conc
org_met_conc = gene_interactions %>%
  ungroup() %>%
  select(Organism, Metabolite, Conc) %>%
  distinct()

# Perform Fisher's Exact Test for each metabolite in each organism
library(foreach)
library(doMC)
registerDoMC(16)

fisher_results = bind_rows(
  foreach(i=1:nrow(org_met_conc)) %dopar% {
  # Select Organism and Metabolite for test
    org_met_conc_test = org_met_conc[i,]
    test_tb = gene_interactions %>%
      filter(
        Organism == org_met_conc_test$Organism,
        Metabolite == org_met_conc_test$Metabolite,
        Conc == org_met_conc_test$Conc
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
      Organism = org_met_conc_test$Organism,
      Metabolite = org_met_conc_test$Metabolite,
      Concentration = org_met_conc_test$Conc,
      Enzyme = frac_enzyme,
      Other = frac_other,
      p = p
    )
  }
)

# Adjust p-values
fisher_results = fisher_results %>%
  mutate(q = p.adjust(p, "BH")) %>%
  arrange(q)

# Save results table
write_tsv(
  fisher_results, "results/Fisher_exact_test_for_enzyme_interactions.tab"
)
