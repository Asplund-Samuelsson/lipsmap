options(width=110)
library(tidyverse)

# Read input file and output directory
args = commandArgs(trailingOnly=TRUE)
infile = args[1]
outdir = args[2]

# Load data
lipsmap = read_tsv(infile) %>%
  # Change significance
  mutate(Sign = ifelse(adj.pvalue < 0.01, "sign", "unsign"))

kegg_modules = read_tsv("data/KEGGgene_module_organism.tab")

kegg_uniprot = read_tsv("data/KEGGgene_uniprot_organism.tab") %>%
  # Clean up the UniProt ID
  mutate(UniProt_entry = str_remove(UniProt, "up:"))

kegg_ec = read_tsv("data/KEGGgene_EC_organism.tab") %>%
  # Clean up the EC
  mutate(EC = str_remove(EC, "ec:"))

uniprot_ec = kegg_uniprot %>%
  inner_join(kegg_ec) %>%
  select(UniProt_entry, EC)

# Determine interactions
interactions = lipsmap %>%
  # Determine metabolite-interacting proteins; at least one significant peptide
  group_by(Organism, Metabolite, Conc, UniProt_entry) %>%
  summarise(Interaction = "sign" %in% Sign) %>%
  # Clarify concentration description
  mutate(Conc = ifelse(Conc == 1, "Low", "High"))

# Determine GO annotations
go_annotations = lipsmap %>%
  select(UniProt_entry, Gene_ontology_GO) %>%
  distinct() %>%
  # Split GO data into rows
  ungroup() %>%
  mutate(GO = str_split(Gene_ontology_GO, "; ")) %>%
  unnest(cols="GO") %>%
  select(-Gene_ontology_GO) %>%
  # Remove genes that have no GO annotation
  filter(!is.na(GO)) %>%
  # Rename to Annotation
  rename(Annotation = GO)

# Determine EC annotations
ec_annotations = uniprot_ec %>%
  # Rename to Annotation
  rename(Annotation = EC)

# Determine pathway annotations
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
  ) %>% mutate(Pathway = str_remove(Pathway, "\\.$")) %>%
  # Rename to Annotation
  rename(Annotation = Pathway)

# Determine KEGG modules annotations
module_annotations = kegg_uniprot %>%
  # Add modules
  inner_join(kegg_modules) %>%
  # Clean up the Module ID
  separate(Module, c(NA, "Annotation"), sep="_") %>%
  select(UniProt_entry, Annotation) %>%
  distinct()

# Function to perform Fisher's exact test by a provided grouping variable
do_fisher = function(annotations) {

  # Determine pathway interactions
  annotated_interactions = inner_join(interactions, annotations)

  # Determine combinations of Metabolite, Conc, and Annotation
  combinations = annotated_interactions %>%
    ungroup() %>%
    select(Conc, Metabolite, Annotation) %>%
    distinct()

  interaction_counts = inner_join(
    # Calculate number of interactions per Annotation
    annotated_interactions %>%
      group_by(Metabolite, Conc, Annotation) %>%
      summarize(
        Annotation_int = sum(Interaction),
        No_Annotation_int = sum(!Interaction)
      ),
    # Calculate number of proteins and interactions in total per group
    annotated_interactions %>%
      select(Metabolite, Conc, UniProt_entry, Interaction) %>%
      distinct() %>%
      group_by(Metabolite, Conc) %>%
      summarize(
        Tot_prot = length(UniProt_entry),
        Tot_int = sum(Interaction)
      )
  ) %>%
  inner_join(
    # Calculate number of proteins with each Annotation
    annotated_interactions %>%
      select(Metabolite, Conc, UniProt_entry, Annotation) %>%
      group_by(Metabolite, Conc, Annotation) %>%
      summarise(Has_Annotation = length(UniProt_entry))
  ) %>%
  mutate(
    # Calculate values needed for contingency table in Fisher's exact test
    No_Annotation = Tot_prot - Has_Annotation,
    Interaction_Annotation = Annotation_int,
    Interaction_Other = Tot_int - Annotation_int,
    No_interaction_Annotation = No_Annotation_int,
    No_interaction_Other = No_Annotation - Tot_int + Annotation_int
  )

  # Perform Fisher's Exact Test per Organism, Metabolite, Conc, and GO
  library(foreach)
  library(doMC)
  registerDoMC(16)

  fisher_results = bind_rows(
    foreach(i=1:nrow(combinations)) %dopar% {
      # Select data for test
      test_tb = inner_join(
        interaction_counts, combinations[i,],
        # Specify "by" to suppress messages
        by = c("Metabolite", "Conc", "Annotation")
      )
      # Set up contingency table
      c_table = matrix(
        c(
          test_tb$Interaction_Annotation, test_tb$Interaction_Other,
          test_tb$No_interaction_Annotation, test_tb$No_interaction_Other
        ), ncol = 2
      )
      # Calculate interaction fraction
      frac_annotation = test_tb$Interaction_Annotation /
        (test_tb$Interaction_Annotation + test_tb$No_interaction_Annotation)
      frac_other = test_tb$Interaction_Other /
        (test_tb$Interaction_Other + test_tb$No_interaction_Other)
      # Perform test
      p = fisher.test(c_table)$p.value
      # Return results
      tibble(
        Metabolite = test_tb$Metabolite,
        Concentration = test_tb$Conc,
        Annotation = test_tb$Annotation,
        Annotation_frac = frac_annotation,
        Other_frac = frac_other,
        p = p
      )
    }
  )

  # Adjust p-values
  fisher_results %>% mutate(q = p.adjust(p, "BH")) %>% arrange(q)

}

# Analyze GO
go_fisher = do_fisher(go_annotations)
write_tsv(
  go_fisher, file.path(outdir, "Fisher_exact_test_for_GO_interactions.tab")
)

# Analyze EC
ec_fisher = do_fisher(ec_annotations)
write_tsv(
  ec_fisher, file.path(outdir, "Fisher_exact_test_for_EC_interactions.tab")
)

# Analyze Pathway
pw_fisher= do_fisher(pathway_annotations)
write_tsv(
  pw_fisher, file.path(outdir, "Fisher_exact_test_for_pathway_interactions.tab")
)

# Analyze Module
mo_fisher = do_fisher(module_annotations)
write_tsv(
  mo_fisher, file.path(outdir, "Fisher_exact_test_for_module_interactions.tab")
)
