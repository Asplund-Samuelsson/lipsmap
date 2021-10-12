options(width=110)
library(ggpubr)
library(tidyverse)

# Load data
lipsmap = read_tsv("data/annotated_comparison_results.tab.gz") %>%
  # Change significance
  mutate(Sign = ifelse(adj.pvalue < 0.01, "sign", "unsign"))

kegg_ec = read_tsv("data/KEGGgene_EC_organism.tab")

kegg_uniprot = read_tsv("data/KEGGgene_uniprot_organism.tab") %>%
  # Clean up the UniProt ID
  mutate(UniProt_entry = str_remove(UniProt, "up:"))

org_uniprot = read_tsv(
  "data/organism_uniprot.tab",
  col_names = c("Organism", "UniProt_entry")
)

cbb_ec = read_tsv("data/cbb_enzymes.tab")

cbb_drains = read_tsv("data/cbb_drains.tab")

# Determine interactions
interactions = lipsmap %>%
  # Determine metabolite-interacting proteins; at least one significant peptide
  group_by(Organism, Metabolite, Conc, UniProt_entry) %>%
  summarise(Interaction = "sign" %in% Sign) %>%
  # Clarify concentration description
  mutate(Conc = ifelse(Conc == 1, "Low", "High"))

# Add enzymes
plot_data = interactions %>%
  # Add EC
  inner_join(kegg_uniprot) %>%
  inner_join(kegg_ec) %>%
  mutate(EC = str_remove(EC, "ec:")) %>%
  # Add and reduce to CBB and Drains
  inner_join(
    bind_rows(
      cbb_ec %>% mutate(Pathway = "CBB"),
      cbb_drains %>% mutate(Pathway = "Drains")
    )
  ) %>%
  # Remove superfluous columns
  select(-UniProt_entry, -UniProt, -KEGGgene, -EC, -Name) %>%
  # Summarise interaction per enzyme
  group_by(
    Organism, Metabolite, Conc, Abbreviation, Pathway
  ) %>%
  summarise(Interaction = T %in% Interaction)

# Determine measured data
measured_data = interactions %>%
  ungroup() %>%
  select(Organism, Metabolite, Conc) %>%
  distinct()

# Determine existing categories
existing_enzymes = org_uniprot %>%
  # Add EC
  inner_join(kegg_uniprot) %>%
  inner_join(kegg_ec) %>%
  mutate(EC = str_remove(EC, "ec:")) %>%
  # Add and reduce to CBB and Drains
  inner_join(bind_rows(cbb_ec, cbb_drains)) %>%
  select(Organism, Abbreviation) %>%
  distinct()

# Determine all possible data to be plotted
possible_data = expand_grid(
  Organism = unique(lipsmap$Organism),
  Abbreviation = unique(plot_data$Abbreviation),
  Conc = c("High", "Low"),
  Metabolite = unique(lipsmap$Metabolite)
)

# Classify missing data
missing_data = possible_data %>%
  anti_join(plot_data) %>%
  # Determine if Enzyme exists in Organism
  left_join(existing_enzymes %>% mutate(Exists = T)) %>%
  mutate(Exists = ifelse(is.na(Exists), F, T)) %>%
  # Determine if Metabolite was measured
  left_join(measured_data %>% mutate(Measured = T)) %>%
  mutate(Measured = ifelse(is.na(Measured), F, T)) %>%
  # Determine Data status
  mutate(
    Data = case_when(
      !Exists ~ "Inexistent",
      !Measured ~ "Unmeasured",
      T ~ "Undetected"
    )
  ) %>%
  # Add back Pathway
  inner_join(
    plot_data %>%
      ungroup() %>%
      select(Abbreviation, Pathway) %>%
      distinct()
  ) %>%
  # Clean up
  select(Organism, Metabolite, Conc, Abbreviation, Pathway, Data) %>%
  mutate(Interaction = F)

# Order organisms
organisms = c("Hydrogenophaga", "Cupriavidus", "Synechococcus", "Synechocystis")
organcols = c("#762a83", "#9970ab","#5aae61","#1b7837")

# Prepare heatmap table
interactions_heat = plot_data %>%
  bind_rows(missing_data) %>%
  mutate(
    Data = factor(
      ifelse(is.na(Data), "", Data),
      levels = c("", "Undetected", "Unmeasured", "Inexistent")
    ),
    Organism = factor(Organism, levels = organisms),
    Abbreviation = factor(
      Abbreviation,
      levels = plot_data %>%
        group_by(Abbreviation) %>%
        summarise(Interaction = sum(Interaction)) %>%
        arrange(-Interaction) %>%
        pull(Abbreviation) %>%
        unique()
    ),
    Metabolite = factor(
      Metabolite,
      levels = plot_data %>%
        group_by(Metabolite) %>%
        summarise(Interaction = sum(Interaction)) %>%
        arrange(-Interaction) %>%
        pull(Metabolite) %>%
        unique()
    ),
    Metabolism = factor(
      ifelse(str_starts(Organism, "Syn"), "Photoautotroph", "Lithoautotroph"),
      levels=c("Photoautotroph", "Lithoautotroph")
    ),
    Friendship = factor(
      ifelse(Organism %in% c("Cupriavidus", "Synechocystis"), "Old", "New"),
      levels=c("Old", "New")
    )
  )

# Plot heatmap
make_heatmap = function(interactions_heat_sub, ytext=T, xtext=T){
  gp = ggplot(
    interactions_heat_sub,
    aes(x=Friendship, y=Metabolism, fill=Organism, alpha=Interaction, shape=Data)
  )
  gp = gp + geom_tile()
  gp = gp + geom_point(mapping = aes(color=Organism), alpha=0.5)
  gp = gp + facet_grid(Abbreviation~Metabolite, switch="both")
  gp = gp + theme_bw()
  gp = gp + scale_fill_manual(values=organcols)
  gp = gp + scale_color_manual(values=organcols)
  gp = gp + scale_alpha_manual(values=c(0,1))
  gp = gp + scale_shape_manual(values=c(NA, 1, 4, 0), drop=F)
  gp = gp + theme(
    strip.background = element_blank(),
    aspect.ratio = 1,
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    strip.text.y.left = element_text(angle=0, hjust=1, vjust=0.5),
    strip.text.x.bottom = element_text(angle=90, hjust=1, vjust=1),
    panel.spacing = unit(0.1, "lines"),
    panel.border = element_rect(color="#d9d9d9")
  )
  if(!ytext){gp = gp + theme(strip.text.y.left = element_blank())}
  if(!xtext){gp = gp + theme(strip.text.x.bottom = element_blank())}
  return(gp)
}

outfile = "results/cbb_drains.pdf"
pdf(outfile, width=11.5, height=6.5, onefile=FALSE)
print(ggarrange(
  make_heatmap(
    interactions_heat %>% filter(Pathway == "CBB", Conc == "High"), T, F
  ),
  make_heatmap(
    interactions_heat %>% filter(Pathway == "CBB", Conc == "Low"), F, F
  ),
  make_heatmap(
    interactions_heat %>% filter(Pathway == "Drains", Conc == "High"), T, T
  ),
  make_heatmap(
    interactions_heat %>% filter(Pathway == "Drains", Conc == "Low"), F, T
  ),
  common.legend=T,
  nrow = 2,
  ncol = 2,
  widths = c(1,1),
  heights = c(1,1),
  align = "hv"
))
dev.off()
