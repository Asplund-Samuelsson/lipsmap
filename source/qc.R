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

# Check how many low concentration interactions are present also in high
interaction = lipsmap %>%
  # Clarify concentration
  mutate(Conc = ifelse(Conc == 1, "Low", "High")) %>%
  # Determine interacting proteins
  group_by(Organism, Metabolite, Conc, UniProt_entry) %>%
  summarise(Interaction = "sign" %in% Sign) %>%
  # Get column with Low and High interaction
  spread(Conc, Interaction) %>%
  # Keep only low interacting proteins
  filter(Low) %>%
  # Count proteins
  group_by(Organism, Metabolite, High) %>%
  summarise(Proteins = length(UniProt_entry)) %>%
  rename(Persists = High) %>%
  mutate(Fraction = Proteins / sum(Proteins))

# Generate percent labels
library(scales)

percent_labels = interaction %>%
  ungroup() %>%
  complete(
    expand(., Organism, Metabolite, Persists), fill=list(Proteins=0, Fraction=0)
  ) %>%
  group_by(Organism, Metabolite) %>%
  filter(sum(Proteins) > 0) %>%
  mutate(
    Percent = paste(
      percent(Fraction), " of ", sum(Proteins), sep=""
    ),
    Fraction = 0.02
  ) %>%
  filter(Persists)

# Organism colours
organisms = rev(
  c("Hydrogenophaga", "Cupriavidus", "Synechococcus", "Synechocystis")
)
organcols = rev(c("#762a83", "#9970ab","#5aae61","#1b7837"))
names(organcols) = organisms

# Order Organism and Metabolite
metabolites = interaction %>%
  group_by(Metabolite) %>%
  summarise(Proteins = sum(Proteins)) %>%
  arrange(-Proteins) %>%
  pull(Metabolite)

interaction = interaction %>%
  mutate(
    Organism = factor(Organism, levels=organisms),
    Metabolite = factor(Metabolite, levels=metabolites)
  )

percent_labels = percent_labels %>%
  mutate(
    Organism = factor(Organism, levels=organisms),
    Metabolite = factor(Metabolite, levels=metabolites)
  )

# Plot it
gp = ggplot(
  interaction,
  aes(y=Fraction, x=Metabolite, fill=Organism, alpha=Persists)
)
gp = gp + geom_col(position="stack")
gp = gp + geom_text(
  data=percent_labels, aes(label=Percent), color="white", size=2,
  hjust=0
)
gp = gp + scale_fill_manual(values=organcols, guide="none")
gp = gp + scale_alpha_manual(values=c(0.5,1), guide="none")
gp = gp + theme_bw()
gp = gp + facet_grid(. ~ Organism)
gp = gp + theme(
  axis.ticks = element_line(colour="black"),
  axis.text = element_text(colour="black"),
#  axis.text.y = element_text(angle=90, hjust=1, vjust=0.5),
  axis.ticks.x = element_blank(),
  axis.text.x = element_blank(),
  strip.background = element_blank(),
  legend.position="top"
)
gp = gp + ylab(
  "Persistence of low concentration interactions in high concentration"
)
gp = gp + coord_flip()

ggsave(file.path(outdir, "persistence.pdf"), gp, w=5.35, h=4)
