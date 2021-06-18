library(tidyverse)
options(width=110)

# Load data
kegg_uniprot = read_tsv("data/KEGGgene_uniprot_organism.tab")
kegg_modules = read_tsv("data/KEGGgene_module_organism.tab")
orthologs = read_tsv(
  "data/uniprot_eggNOG.tab",
  col_names = c("UniProt_entry", "Ortholog")
)
eggnog_category_descriptions = read_tsv(
  "data/ortholog_category_descriptions.tab"
)
module_description = read_tsv(
  "data/module_description.tab", col_names=c("Module", "Description")
)
eggnog_annotations = read_tsv("data/eggNOG_annotations.tab")

# Define organisms and colors
organisms = c("Hydrogenophaga", "Cupriavidus", "Synechococcus", "Synechocystis")
organcols = c("#762a83", "#9970ab","#5aae61","#1b7837")

# Determine eggNOG annotations
eggnog_annotations_unique = eggnog_annotations %>%
  # Select relevant data
  select(Ortholog, Version, Category) %>%
  # Split multiple categories between each character
  mutate(Category = str_split(Category, "")) %>%
  # Unnest categories per ortholog into rows
  unnest(cols="Category") %>%
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

# Combine
module_category = kegg_uniprot %>%
  inner_join(kegg_modules) %>%
  mutate(UniProt_entry = str_remove(UniProt, "up:")) %>%
  separate(Module, c(NA, "Module"), "_") %>%
  inner_join(orthologs) %>%
  inner_join(eggnog_annotations_unique) %>%
  inner_join(eggnog_category_descriptions) %>%
  mutate(Label = paste("[", Category, "] ", Short_description, sep="")) %>%
  rename(Ortholog_description = Label) %>%
  select(-Description) %>%
  inner_join(
    module_description %>% mutate(Module = str_remove(Module, "md:"))
  ) %>%
  rename(Module_description = Description) %>%
  group_by(Ortholog_description, Module_description, Organism) %>%
  summarise(Orthologs = length(Ortholog)) %>%
  ungroup() %>%
  complete(
    Ortholog_description, Module_description, Organism,
    fill = list(Orthologs = 0)
  ) %>%
  mutate(Organism = factor(Organism, levels=organisms))

# Filter to top 50 modules
wanted_modules = module_category %>%
  group_by(Module_description) %>%
  summarise(Orthologs = sum(Orthologs)) %>%
  top_n(50, Orthologs) %>%
  arrange(-Orthologs) %>%
  pull(Module_description)

module_category = module_category %>%
  filter(Module_description %in% wanted_modules) %>%
  mutate(Module_description = factor(Module_description, levels=wanted_modules))

# Plot it
gp = ggplot(
  module_category,
  aes(x=Organism, y=Orthologs, group=Organism, fill=Organism)
)
gp = gp + geom_col(position=position_dodge())
gp = gp + facet_grid(Module_description~Ortholog_description, scales="free_y")
gp = gp + scale_fill_manual(values=organcols)
gp = gp + theme_bw()
gp = gp + theme(
  strip.background = element_blank(),
  axis.text = element_text(color="black"),
  axis.ticks = element_line(color="black"),
  strip.text.y = element_text(angle=0, hjust=0, vjust=0.5),
  strip.text.x = element_text(angle=90, hjust=0, vjust=0.5),
  axis.text.x = element_blank(),
  axis.title.x = element_blank(),
  axis.ticks.x = element_blank(),
  panel.grid = element_blank()
)

ggsave("results/category_module_overlap.pdf", gp, w=18, h=28)
