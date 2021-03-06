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

orthologs = read_tsv(
  "data/uniprot_eggNOG.tab",
  col_names = c("UniProt_entry", "Ortholog")
)

org_uniprot = read_tsv(
  "data/organism_uniprot.tab",
  col_names = c("Organism", "UniProt_entry")
)

uniprot_locus = read_tsv(
    "data/uniprot_locus_complete.tab", col_names=c("UniProt_entry", "Locus")
)

eggnog_annotations = read_tsv("data/eggNOG_annotations.tab")

eggnog_category_descriptions = read_tsv(
  "data/ortholog_category_descriptions.tab"
)

kegg_uniprot = read_tsv("data/KEGGgene_uniprot_organism.tab") %>%
  # Clean up the UniProt ID
  mutate(UniProt_entry = str_remove(UniProt, "up:"))

kegg_modules = read_tsv("data/KEGGgene_module_organism.tab")

kegg_ec = read_tsv("data/KEGGgene_EC_organism.tab")

modules = kegg_uniprot %>%
  # Add modules
  inner_join(kegg_modules) %>%
  # Clean up the Module ID
  separate(Module, c(NA, "Module"), sep="_") %>%
  select(UniProt_entry, Module) %>%
  distinct()

module_description = read_tsv(
  "data/module_description.tab", col_names=c("Module", "Description")
)

module_description = module_description %>%
  mutate(Module = str_remove(Module, "md:"))

module_description_short = read_tsv("data/module_description_short.tab") %>%
  mutate(Label = paste(Short_description, " [", Module, "]", sep=""))

ec_description = read_tsv(
  "data/EC_description.tab", col_names=c("EC", "Description")
)

# Correct Ortholog names based on eggNOG version
eggnog_annotations = eggnog_annotations %>%
  # Add missing version string to non-COG orthologs
  mutate(
    Ortholog = ifelse(
      str_starts(Ortholog, "COG"),
      Ortholog,
      ifelse(
        Version == 4.1,
        paste("ENOG41", Ortholog, sep=""),
        paste("ENOG50", Ortholog, sep="")
      )
    )
  ) %>%
  # Remove version and make rows unique
  select(-Version) %>%
  distinct()

eggnog_description = eggnog_annotations %>%
  select(Ortholog, Description) %>%
  distinct()

# Determine eggNOG annotations
eggnog_annotations_unique = eggnog_annotations %>%
  # Select relevant data
  select(Ortholog, Category) %>%
  # Split multiple categories between each character
  mutate(Category = str_split(Category, "")) %>%
  # Unnest categories per ortholog into rows
  unnest(cols="Category") %>%
  # Make rows unique
  distinct()

# Determine interactions
interactions = lipsmap %>%
  # Determine metabolite-interacting proteins; at least one significant peptide
  group_by(Organism, Metabolite, Conc, UniProt_entry) %>%
  summarise(Interaction = "sign" %in% Sign) %>%
  # Clarify concentration description
  mutate(Conc = ifelse(Conc == 1, "Low", "High"))

# Determine experiments
experiments = interactions %>%
  select(Organism, Metabolite, Conc) %>%
  distinct() %>%
  # Expand to all proteins per Organism
  inner_join(org_uniprot)

# Add missing data to interactions
interactions = left_join(experiments, interactions)

# Add annotation data
annotated_interactions = interactions %>%
  ungroup() %>%
  # Add KEGG gene names
  left_join(kegg_uniprot) %>%
  select(-UniProt) %>%
  # Add KEGG modules and description
  left_join(modules) %>%
  left_join(module_description) %>%
  rename(Module_description = Description) %>%
  # Add KEGG EC and description
  left_join(kegg_ec) %>%
  left_join(ec_description) %>%
  rename(EC_description = Description) %>%
  # Add eggNOG orthologs, category, and description
  left_join(orthologs) %>%
  left_join(eggnog_annotations_unique) %>%
  left_join(eggnog_category_descriptions) %>%
  rename(
    Ortholog_category = Description,
    Ortholog_short_category = Short_description
  ) %>%
  left_join(
    eggnog_description %>%
      group_by(Ortholog) %>%
      summarise(Ortholog_description = paste(Description, collapse="; "))
  ) %>%
  # Add locus ID
  left_join(uniprot_locus) %>%
  # Clean up
  select(-KEGGgene) %>%
  mutate(
    EC = str_remove(EC, "ec:"),
    EC_description = unlist(lapply(str_split(EC_description, ";"), "[[", 1))
  ) %>%
  select(
    Organism, Metabolite, Conc, UniProt_entry, Locus, Interaction,
    EC, EC_description, Module, Module_description,
    Ortholog, Ortholog_description, Ortholog_category
  )

# Save Excel file
library(openxlsx)

wb = createWorkbook()
modifyBaseFont(wb, fontSize = 11, fontColour = "black", fontName = "Arial")
addWorksheet(wb, "Ortholog EC module interactions")
writeData(
  wb, "Ortholog EC module interactions",
  annotated_interactions,
  headerStyle = createStyle(textDecoration = "bold"),
  keepNA=T, na.string="NA"
)
setColWidths(
  wb, "Ortholog EC module interactions",
  cols=1:ncol(annotated_interactions), widths = "auto"
)
freezePane(wb, "Ortholog EC module interactions", firstRow = TRUE)
saveWorkbook(
  wb, file=file.path(outdir, "ortholog_ec_module_interactions.xlsx"),
  overwrite=T
)
