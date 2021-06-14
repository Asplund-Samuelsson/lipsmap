options(width=110)
library(tidyverse)

# Load data
lipsmap = read_tsv(
  "data/annotated_comparison_results.tab.gz",
  col_types = cols(Protein.names = col_character())
)

kegg_uniprot = read_tsv("data/KEGGgene_uniprot_organism.tab")

kegg_modules = read_tsv("data/KEGGgene_module_organism.tab")

module_reaction = read_tsv(
  "data/module_reaction.tab", col_names=c("Module", "Reaction")
)

reaction_compound = read_tsv(
  "data/reaction_compound.tab", col_names=c("Reaction", "Compound")
)

metabolite_compound = read_tsv(
  "data/compound_kegg.tab", col_names=c("Metabolite", "Compound")
)

module_description = read_tsv(
  "data/module_description.tab", col_names=c("Module", "Description")
)

# Define organisms and colors
organisms = c("Hydrogenophaga", "Cupriavidus", "Synechococcus", "Synechocystis")
organcols = c("#762a83", "#9970ab","#5aae61","#1b7837")

# Clean up strings
module_compound = module_reaction %>%
  inner_join(reaction_compound) %>%
  select(-Reaction) %>%
  mutate(
    Module = str_remove(Module, "md:"),
    Compound = str_remove(Compound, "cpd:")
  ) %>%
  distinct()

module_description = module_description %>%
  mutate(Module = str_remove(Module, "md:"))

modules = kegg_uniprot %>%
  # Clean up the UniProt ID
  mutate(UniProt_entry = str_remove(UniProt, "up:")) %>%
  # Add modules
  inner_join(kegg_modules) %>%
  # Clean up the Module ID
  separate(Module, c(NA, "Module"), sep="_") %>%
  select(UniProt_entry, Module) %>%
  distinct()

# Determine interactions
interactions = lipsmap %>%
  # Determine metabolite-interacting proteins; at least one significant peptide
  group_by(Organism, Metabolite, Conc, UniProt_entry) %>%
  summarise(Interaction = "sign" %in% Sign) %>%
  # Clarify concentration description
  mutate(Conc = ifelse(Conc == 1, "Low", "High"))

interactions = interactions %>%
  # Add Module to Interaction
  inner_join(modules) %>%
  # Add Compound to Metabolite
  inner_join(metabolite_compound) %>%
  # Determine if Compound is in Module
  left_join(module_compound %>% mutate(Relation = "Internal")) %>%
  mutate(Relation = replace_na(Relation, "External"))

# Count Interaction per Module
interaction_counts = interactions %>%
  group_by(Organism, Metabolite, Conc, Module, Relation) %>%
  summarise(Interactions = sum(Interaction))

# Add Description and save output
write_tsv(
  interaction_counts %>%
    arrange(-Interactions) %>%
    inner_join(module_description),
  "results/module_interaction_summary.tab"
)

# Add missing data
interaction_counts = interaction_counts %>%
  ungroup() %>%
  complete(
    Organism, Metabolite, Conc, Module, fill = list(Interactions = 0)
  ) %>%
  # Add the KEGG Compound IDs again
  inner_join(metabolite_compound) %>%
  # Determine if Compound is in Module
  select(-Relation) %>%
  left_join(module_compound %>% mutate(Relation = "Internal")) %>%
  mutate(Relation = replace_na(Relation, "External"))

# Order variables
interaction_counts = interaction_counts %>%
  inner_join(module_description) %>%
  mutate(
    Metabolite = factor(
      Metabolite,
      levels = (
        group_by(., Metabolite) %>%
          summarise(Interactions = sum(Interactions)) %>%
          arrange(-Interactions) %>%
          pull(Metabolite)
      )
    ),
    Description = factor(
      Description,
      levels = (
        group_by(., Description) %>%
          summarise(Interactions = sum(Interactions)) %>%
          arrange(-Interactions) %>%
          pull(Description)
      )
    )
  )

interaction_counts = mutate(
  interaction_counts, Organism = factor(Organism, levels = organisms)
)

# Make plot
gp = ggplot(
  # Consider only top modules
  filter(interaction_counts, as.numeric(Description) %in% 1:40),
  aes(
    x=Organism, alpha=Conc,
    group=factor(
      paste(Organism, Conc),
      levels=paste(rep(organisms, each=2), rep(c("High", "Low"), 4))
    ),
    fill=Organism, y=Interactions
  )
)
gp = gp + geom_rect(
    aes(linetype=Relation),
    xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, colour="black", size=0.4,
    alpha=0.6
  )
gp = gp + geom_col(position=position_dodge(width=1))
gp = gp + theme_bw()
gp = gp + facet_grid(Description~Metabolite)
gp = gp + theme(
  strip.background = element_blank(),
  strip.text.y = element_text(angle=0, hjust=0, vjust=0.5),
  strip.text.x = element_text(angle=90, hjust=0, vjust=0.5),
  axis.ticks = element_blank(),
  axis.text = element_blank(),
  panel.grid = element_blank(),
  aspect.ratio=1,
  panel.border=element_blank(),
  axis.title.x=element_blank(),
  legend.position="top"
)
gp = gp + scale_alpha_manual(values=c(1,0.25))
gp = gp + scale_fill_manual(values=organcols)
gp = gp + scale_linetype_manual(values=c(3,1))

ggsave("results/module_interactions.pdf", gp, w=15, h=16)
