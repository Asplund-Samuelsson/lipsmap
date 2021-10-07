options(width=110)
library(tidyverse)
library(ggarrange)

# Load data
lipsmap = read_tsv("data/annotated_comparison_results.tab.gz") %>%
  # Change significance
  mutate(Sign = ifelse(adj.pvalue < 0.01, "sign", "unsign"))

ccm = read_csv("data/CCM_regulatory_proteins.csv")

# Clean up CCM annotations
ccm = ccm %>%
  select(Gene_name, Locus_tag, Category, Protein_complex_1) %>%
  bind_rows(
    ccm %>%
      select(Gene_name, Locus_tag, Category, Protein_complex_2) %>%
      filter(!is.na(Protein_complex_2)) %>%
      rename(Protein_complex_1 = Protein_complex_2)
  ) %>%
  rename(Name = Gene_name, Locus = Locus_tag, Complex = Protein_complex_1)

# Determine interactions
interactions = lipsmap %>%
  # Select only Synechocystis
  filter(Organism == "Synechocystis") %>%
  # Split Locus
  select(Metabolite, Conc, Locus, Sign) %>%
  mutate(Locus = str_split(str_remove(Locus, ";"), " ")) %>%
  unnest(Locus) %>%
  # Determine metabolite-interacting proteins; at least one significant peptide
  group_by(Metabolite, Conc, Locus) %>%
  summarise(Interaction = "sign" %in% Sign) %>%
  # Clarify concentration description
  mutate(Conc = ifelse(Conc == 1, "Low", "High"))

# Determine CCM interactions
ccm_interactions = inner_join(ccm, interactions)

# Determine all possible data points
possible_data = expand_grid(
  Locus = unique(ccm$Locus),
  Metabolite = unique(interactions$Metabolite),
  Conc = c("High", "Low")
)

# Determine missing CCM proteins
missing_data = possible_data %>%
  anti_join(ccm_interactions) %>%
  inner_join(ccm) %>%
  mutate(Interaction = F, Data = "Undetected")

# Determine order of wanted metabolites
metabolites = ccm_interactions %>%
  group_by(Metabolite) %>%
  summarise(Interaction = sum(Interaction)) %>%
  arrange(-Interaction) %>%
  filter(Interaction > 0) %>%
  pull(Metabolite)

# Prepare plotting data
plot_data = bind_rows(
  ccm_interactions %>% mutate(Data = ""),
  missing_data
) %>%
  # Select only wanted metabolites
  filter(Metabolite %in% metabolites) %>%
  mutate(
    # Prepare a Gene label
    Gene = paste(Name, " [", Locus, "]", sep=""),
    # Order concentrations
    Conc = factor(Conc, levels=c("High", "Low")),
    # Order metabolites
    Metabolite = factor(Metabolite, levels=metabolites),
    # Replace missing Complex by creating new Category
    Category = ifelse(
      is.na(Complex),
      Category,
      paste(Category, " (", Complex, ")", sep="")
    )
  ) %>%
  select(-Name, -Locus, -Complex)

plot_data = plot_data %>%
  # Order genes
  mutate(
    Gene = factor(
      Gene,
      plot_data %>%
        group_by(Gene) %>%
        summarise(Interaction = sum(Interaction)) %>%
        arrange(Interaction) %>%
        pull(Gene)
    )
  )

# Make plot
gp = ggplot(plot_data, aes(x=Metabolite, y=Gene, shape=Data))
gp = gp + geom_tile(fill="#1b7837", mapping=aes(alpha=Interaction))
gp = gp + geom_point(alpha=0.5, color="#1b7837")
gp = gp + scale_shape_manual(values=c(NA, 1))
gp = gp + scale_alpha_manual(values=c(0,1))
gp = gp + facet_grid(Category~Conc, scales="free_y", space="free")
gp = gp + theme_bw()
gp = gp + theme(
  strip.background = element_blank(),
  axis.ticks = element_line(color="black"),
  axis.text = element_text(color="black"),
  axis.text.x = element_text(angle=90, hjust=1, vjust=1),
  legend.position = "top",
  legend.title = element_blank()
)
gp = gp + guides(alpha=F)

ggsave("results/ccm.pdf", gp, w=8, h=8)
