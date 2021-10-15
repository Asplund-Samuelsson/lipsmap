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

# Determine interactions
interactions = lipsmap %>%
  # Determine metabolite-interacting proteins; at least one significant peptide
  group_by(Organism, Metabolite, Conc, UniProt_entry) %>%
  mutate(Interaction = "sign" %in% Sign) %>%
  # Clarify concentration description
  mutate(Conc = ifelse(Conc == 1, "Low", "High")) %>%
  # Select relvant data
  select(Organism, Metabolite, Conc, Peptide_ID, UniProt_entry, Interaction)

# Order organisms
organisms = rev(
  c("Hydrogenophaga", "Cupriavidus", "Synechococcus", "Synechocystis")
)

# Summarise number of peptides per protein
peptides = interactions %>%
  group_by(Organism, Metabolite, Conc, UniProt_entry, Interaction) %>%
  summarise(Peptides = length(Peptide_ID)) %>%
  mutate(Organism = factor(Organism, levels=organisms))

# Add missing data
possible_data = expand_grid(
  Organism = unique(peptides$Organism),
  Metabolite = unique(peptides$Metabolite),
  Conc = unique(peptides$Conc),
  Interaction = unique(peptides$Interaction)
)

missing_data = possible_data %>%
  anti_join(peptides) %>%
  mutate(UniProt_entry = NA, Peptides = 0)

plot_data = bind_rows(peptides, missing_data)

# Determine order of metabolites
metabolites = ungroup(plot_data) %>%
  group_by(Metabolite) %>%
  summarise(Peptides = sum(ifelse(Interaction, Peptides, 0))) %>%
  arrange(-Peptides) %>%
  pull(Metabolite)

plot_data = plot_data %>%
  # Sort Metabolite by number of Interactions
  ungroup() %>%
  mutate(
    Metabolite = factor(Metabolite, levels = metabolites),
    Conc = ifelse(Conc == "High", "High concentration", "Low concentration")
  )

# Plot it
gp = ggplot(
  plot_data,
  aes(
    x=Metabolite, y=Peptides,
    fill=Interaction, colour=Interaction,
    group=paste(Metabolite, Interaction)
  )
)
gp = gp + geom_boxplot(outlier.size=0.1, outlier.stroke=0.1)
gp = gp + facet_grid(Organism ~ Conc)
gp = gp + scale_fill_manual(values=c("#e08214", "#8073ac"))
gp = gp + scale_color_manual(values=c("#7f3b08", "#2d004b"))
gp = gp + theme_bw()
gp = gp + theme(
  axis.ticks = element_line(color="black"),
  axis.text = element_text(color="black"),
  strip.background = element_blank(),
  axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
  legend.position = "top"
)
gp = gp + ylab("Peptides per protein")

ggsave(file.path(outdir,"peptides_per_protein.pdf"), gp, w=9, h=11)


# Count number of detected peptides per experiment
tot_peptides = peptides %>%
  group_by(Organism, Metabolite, Conc) %>%
  summarise(Peptides = sum(Peptides))

# Count number of significant peptides per experiment
tot_significant = tot_peptides %>%
  select(-Peptides) %>%
  left_join(
    lipsmap %>%
      group_by(Organism, Metabolite, Conc) %>%
      summarise(Peptides = sum(Sign == "sign")) %>%
      mutate(Conc = ifelse(Conc == 1, "Low", "High"))
  ) %>%
  mutate(Organism = factor(Organism, levels=organisms))

# Organism colours
organcols = rev(c("#762a83", "#9970ab","#5aae61","#1b7837"))
names(organcols) = organisms

# Plot it
gp = ggplot(
  tot_peptides,
  aes(y=Peptides, x=Metabolite, fill=Organism, group=Conc, alpha=Conc)
)
gp = gp + geom_col(position="dodge")
gp = gp + geom_point(
  data=tot_significant, shape=21, alpha=1, position=position_dodge(width=1),
  fill="white", colour="black", size=1
)
gp = gp + scale_fill_manual(values=organcols, guide=F)
gp = gp + scale_alpha_manual(values=c(1,0.5), name="Concentration")
gp = gp + theme_bw()
gp = gp + facet_grid(Organism ~ .)
gp = gp + theme(
  axis.ticks = element_line(colour="black"),
  axis.text = element_text(colour="black"),
  axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
  strip.background = element_blank(),
  legend.position="top"
)
gp = gp + ylab("Number of detected peptides (circles indicate significant)")

ggsave(file.path(outdir, "peptides_per_experiment.pdf"), gp, w=6, h=6)
