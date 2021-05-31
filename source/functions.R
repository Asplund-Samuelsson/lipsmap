options(width=110)
library(tidyverse)

# Load data
lipsmap = read_tsv(
  "data/annotated_comparison_results.tab.gz",
  col_types = cols(Protein.names = col_character())
)

# Determine function interactions
function_interactions = lipsmap %>%
  # Determine metabolite-interacting proteins; at least one significant peptide
  group_by(Organism, Metabolite, Conc, UniProt_entry) %>%
  summarise(Interaction = "sign" %in% Sign) %>%
  # Clarify concentration description
  mutate(Conc = ifelse(Conc == 1, "Low", "High")) %>%
  # Add back gene ontology function data
  inner_join(distinct(select(lipsmap, UniProt_entry, Gene_ontology_GO))) %>%
  # Split GO data into rows
  ungroup() %>%
  mutate(GO = str_split(Gene_ontology_GO, "; ")) %>%
  unnest(cols="GO") %>%
  select(-Gene_ontology_GO)

# Determine combinations of Organism, Metabolite, Conc, and GO
combinations = function_interactions %>%
  ungroup() %>%
  select(Organism, Conc, Metabolite, GO) %>%
  distinct() %>%
  # Do not test when GO is missing (will still count towards total proteins)
  filter(!is.na(GO))

go_interaction_counts = inner_join(
  # Calculate number of interactions per GO
  function_interactions %>%
    group_by(Organism, Metabolite, Conc, GO) %>%
    summarize(
      GO_int = sum(Interaction),
      No_GO_int = sum(!Interaction)
    ),
  # Calculate number of proteins and interactions in total per group
  function_interactions %>%
    select(Organism, Metabolite, Conc, UniProt_entry, Interaction) %>%
    distinct() %>%
    group_by(Organism, Metabolite, Conc) %>%
    summarize(
      Tot_prot = length(UniProt_entry),
      Tot_int = sum(Interaction)
    )
) %>%
inner_join(
  # Calculate number of proteins with each GO
  function_interactions %>%
    select(Organism, Metabolite, Conc, UniProt_entry, GO) %>%
    group_by(Organism, Metabolite, Conc, GO) %>%
    summarise(Has_GO = length(UniProt_entry))
) %>%
mutate(
  No_GO = Tot_prot - Has_GO,
  Interaction_GO = GO_int,
  Interaction_Other = Tot_int - GO_int,
  No_interaction_GO = No_GO_int,
  No_interaction_Other = No_GO - Tot_int + GO_int
)

# Perform Fisher's Exact Test per Organism, Metabolite, Conc, and GO
library(foreach)
library(doMC)
registerDoMC(16)

fisher_results = bind_rows(
  foreach(i=1:nrow(combinations)) %dopar% {
    # Select data for test
    test_tb = inner_join(
      go_interaction_counts, combinations[i,],
      # Specify "by" to suppress messages
      by = c("Organism", "Metabolite", "Conc", "GO")
    )
    # Set up contingency table
    c_table = matrix(
      c(
        test_tb$Interaction_GO, test_tb$Interaction_Other,
        test_tb$No_interaction_GO, test_tb$No_interaction_Other
      ), ncol = 2
    )
    # Calculate interaction fraction
    frac_go = test_tb$Interaction_GO /
      (test_tb$Interaction_GO + test_tb$No_interaction_GO)
    frac_other = test_tb$Interaction_Other /
      (test_tb$Interaction_Other + test_tb$No_interaction_Other)
    # Perform test
    p = fisher.test(c_table)$p.value
    # Return results
    tibble(
      Organism = test_tb$Organism,
      Metabolite = test_tb$Metabolite,
      Concentration = test_tb$Conc,
      GO = test_tb$GO,
      GO_frac = frac_go,
      Other_frac = frac_other,
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
  fisher_results, "results/Fisher_exact_test_for_GO_interactions.tab"
)


# Determine enzyme interactions
enzyme_interactions = lipsmap %>%
  # Determine metabolite-interacting proteins; at least one significant peptide
  group_by(Metabolite, Conc, UniProt_entry) %>%
  summarise(Interaction = "sign" %in% Sign) %>%
  # Clarify concentration description
  mutate(Conc = ifelse(Conc == 1, "Low", "High")) %>%
  # Add back EC data
  inner_join(distinct(select(lipsmap, UniProt_entry, EC_number))) %>%
  # Split EC data into rows
  ungroup() %>%
  mutate(EC = str_split(EC_number, "; ")) %>%
  unnest(cols="EC") %>%
  select(-EC_number) %>%
  # Remove proteins that are not enzymes
  filter(!is.na(EC))

# Determine combinations of Organism, Metabolite, Conc, and EC
combinations = enzyme_interactions %>%
  ungroup() %>%
  select(Conc, Metabolite, EC) %>%
  distinct()

ec_interaction_counts = inner_join(
  # Calculate number of interactions per EC
  enzyme_interactions %>%
    group_by(Metabolite, Conc, EC) %>%
    summarize(
      EC_int = sum(Interaction),
      No_EC_int = sum(!Interaction)
    ),
  # Calculate number of proteins and interactions in total per group
  enzyme_interactions %>%
    select(Metabolite, Conc, UniProt_entry, Interaction) %>%
    distinct() %>%
    group_by(Metabolite, Conc) %>%
    summarize(
      Tot_prot = length(UniProt_entry),
      Tot_int = sum(Interaction)
    )
) %>%
inner_join(
  # Calculate number of proteins with each EC
  enzyme_interactions %>%
    select(Metabolite, Conc, UniProt_entry, EC) %>%
    group_by(Metabolite, Conc, EC) %>%
    summarise(Has_EC = length(UniProt_entry))
) %>%
mutate(
  No_EC = Tot_prot - Has_EC,
  Interaction_EC = EC_int,
  Interaction_Other = Tot_int - EC_int,
  No_interaction_EC = No_EC_int,
  No_interaction_Other = No_EC - Tot_int + EC_int
)

# Perform Fisher's Exact Test per Organism, Metabolite, Conc, and GO
library(foreach)
library(doMC)
registerDoMC(16)

fisher_results = bind_rows(
  foreach(i=1:nrow(combinations)) %dopar% {
    # Select data for test
    test_tb = inner_join(
      ec_interaction_counts, combinations[i,],
      # Specify "by" to suppress messages
      by = c("Metabolite", "Conc", "EC")
    )
    # Set up contingency table
    c_table = matrix(
      c(
        test_tb$Interaction_EC, test_tb$Interaction_Other,
        test_tb$No_interaction_EC, test_tb$No_interaction_Other
      ), ncol = 2
    )
    # Calculate interaction fraction
    frac_ec = test_tb$Interaction_EC /
      (test_tb$Interaction_EC + test_tb$No_interaction_EC)
    frac_other = test_tb$Interaction_Other /
      (test_tb$Interaction_Other + test_tb$No_interaction_Other)
    # Perform test
    p = fisher.test(c_table)$p.value
    # Return results
    tibble(
      Metabolite = test_tb$Metabolite,
      Concentration = test_tb$Conc,
      EC = test_tb$EC,
      EC_frac = frac_ec,
      Other_frac = frac_other,
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
  fisher_results, "results/Fisher_exact_test_for_EC_interactions.tab"
)


# Extract pathway annotations
pathway_annotations = lipsmap %>%
  # Select relevant data
  select(UniProt_entry, Pathway) %>%
  # Remove proteins without a Pathway annotation
  filter(!is.na(Pathway)) %>%
  distinct() %>%
  # Extract the Pathway names
  mutate(Pathway = str_split(Pathway, "PATHWAY: ")) %>%
  unnest(Pathway) %>%
  filter(Pathway != "") %>%
  mutate(Pathway = sapply(str_split(Pathway, "\\. "), "[[", 1)) %>%
  mutate(Pathway = str_split(Pathway, "; ")) %>%
  mutate(
    Pathway = unlist(map(
      Pathway,
      function(x){if (length(x) == 1) {x[[1]]} else {x[[2]]}}
    ))
  ) %>% mutate(Pathway = str_remove(Pathway, "\\.$"))

# Determine pathway interactions
pathway_interactions = lipsmap %>%
  # Determine metabolite-interacting proteins; at least one significant peptide
  group_by(Metabolite, Conc, UniProt_entry) %>%
  summarise(Interaction = "sign" %in% Sign) %>%
  # Clarify concentration description
  mutate(Conc = ifelse(Conc == 1, "Low", "High")) %>%
  # Add back Pathway data
  inner_join(pathway_annotations)

# Determine combinations of Organism, Metabolite, Conc, and Pathway
combinations = pathway_interactions %>%
  ungroup() %>%
  select(Conc, Metabolite, Pathway) %>%
  distinct()

pathway_interaction_counts = inner_join(
  # Calculate number of interactions per Pathway
  pathway_interactions %>%
    group_by(Metabolite, Conc, Pathway) %>%
    summarize(
      Pathway_int = sum(Interaction),
      No_Pathway_int = sum(!Interaction)
    ),
  # Calculate number of proteins and interactions in total per group
  pathway_interactions %>%
    select(Metabolite, Conc, UniProt_entry, Interaction) %>%
    distinct() %>%
    group_by(Metabolite, Conc) %>%
    summarize(
      Tot_prot = length(UniProt_entry),
      Tot_int = sum(Interaction)
    )
) %>%
inner_join(
  # Calculate number of proteins with each Pathway
  pathway_interactions %>%
    select(Metabolite, Conc, UniProt_entry, Pathway) %>%
    group_by(Metabolite, Conc, Pathway) %>%
    summarise(Has_Pathway = length(UniProt_entry))
) %>%
mutate(
  No_Pathway = Tot_prot - Has_Pathway,
  Interaction_Pathway = Pathway_int,
  Interaction_Other = Tot_int - Pathway_int,
  No_interaction_Pathway = No_Pathway_int,
  No_interaction_Other = No_Pathway - Tot_int + Pathway_int
)

# Perform Fisher's Exact Test per Organism, Metabolite, Conc, and GO
library(foreach)
library(doMC)
registerDoMC(16)

fisher_results = bind_rows(
  foreach(i=1:nrow(combinations)) %dopar% {
    # Select data for test
    test_tb = inner_join(
      pathway_interaction_counts, combinations[i,],
      # Specify "by" to suppress messages
      by = c("Metabolite", "Conc", "Pathway")
    )
    # Set up contingency table
    c_table = matrix(
      c(
        test_tb$Interaction_Pathway, test_tb$Interaction_Other,
        test_tb$No_interaction_Pathway, test_tb$No_interaction_Other
      ), ncol = 2
    )
    # Calculate interaction fraction
    frac_pathway = test_tb$Interaction_Pathway /
      (test_tb$Interaction_Pathway + test_tb$No_interaction_Pathway)
    frac_other = test_tb$Interaction_Other /
      (test_tb$Interaction_Other + test_tb$No_interaction_Other)
    # Perform test
    p = fisher.test(c_table)$p.value
    # Return results
    tibble(
      Metabolite = test_tb$Metabolite,
      Concentration = test_tb$Conc,
      Pathway = test_tb$Pathway,
      Pathway_frac = frac_pathway,
      Other_frac = frac_other,
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
  fisher_results, "results/Fisher_exact_test_for_pathway_interactions.tab"
)
