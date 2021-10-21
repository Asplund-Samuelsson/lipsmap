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

# Order metabolites
metabolites_2 = lipsmap %>%
  group_by(Metabolite) %>%
  summarise(Significant = sum(Sign == "sign", na.rm=T)) %>%
  arrange(-Significant) %>%
  pull(Metabolite)

# Make Hydrogenophaga metabolites come first
metabolites_2 = c(
  metabolites_2[which(metabolites_2 %in% (lipsmap %>% filter(Organism == "Hydrogenophaga") %>% pull(Metabolite) %>% unique()))],
  metabolites_2[which(!(metabolites_2 %in% (lipsmap %>% filter(Organism == "Hydrogenophaga") %>% pull(Metabolite) %>% unique())))]
)

# Investigate fold change and q value of peptides going from low to high
peptides = lipsmap %>%
  select(
    Organism, Metabolite, Conc, log2FC,
    adj.pvalue, UniProt_entry, Peptide,
    Sign
  ) %>%
  mutate(
    Conc = ifelse(Conc == 1, "Low", "High"),
    Significant = (Sign == "sign")
  ) %>%
  rename(q = adj.pvalue) %>%
  select(-Sign) %>%
  distinct() %>%
  mutate(
    Organism = factor(Organism, levels=organisms),
    Metabolite = factor(
      Metabolite,
      levels = metabolites_2
    )
  ) %>%
  filter(!is.na(Significant))

# Make comparison of log2FC in Low and High
peptide_fc = peptides %>%
  select(-q) %>%
  spread(Conc, log2FC) %>%
  filter(!(is.na(High) | is.na(Low)))

# Determine limits
fc_lims = range(c(peptide_fc$High, peptide_fc$Low))

# Correlate peptide log2FC
peptide_fc_cor = peptide_fc %>%
  group_by(Organism, Metabolite, Significant) %>%
  summarise(Correlation = cor(High, Low, method="spearman")) %>%
  filter(!is.na(Correlation)) %>%
  mutate(
    High = ifelse(Significant, fc_lims[1] + 2, fc_lims[2] - 2),
    Low = ifelse(Significant, fc_lims[2] - 1, fc_lims[1] + 1),
    Label = round(Correlation, 3)
  )

# Plot it
make_plot = function (mets, x_title=T, y_title=T) {
  gp = ggplot(
    peptide_fc %>% filter(!Significant, Metabolite %in% mets),
    aes(x=Low, y=High, color=Organism)
  )
  gp = gp + geom_abline(intercept=0, slope=1, size=0.2, color = "grey")
  gp = gp + geom_point(size=0.5, stroke=0)
  gp = gp + geom_point(
    data=peptide_fc %>% filter(Significant, Metabolite %in% mets),
    color="black", size=0.5, stroke=0
  )
  gp = gp + geom_text(
    data=peptide_fc_cor %>% filter(Significant, Metabolite %in% mets),
    mapping=aes(label=Label), color="black", size=2.8,
    hjust=1
  )
  gp = gp + geom_text(
    data=peptide_fc_cor %>% filter(!Significant, Metabolite %in% mets),
    mapping=aes(label=Label), size=2.8,
    hjust=0
  )
  gp = gp + scale_color_manual(values=organcols)
  gp = gp + scale_x_continuous(n.breaks=3, limits=fc_lims)
  gp = gp + scale_y_continuous(n.breaks=3, limits=fc_lims)
  gp = gp + facet_grid(Organism ~ Metabolite)
  gp = gp + theme_bw()
  gp = gp + theme(
    axis.ticks = element_line(color="black"),
    axis.text = element_text(color="black", size=8),
    strip.background = element_blank(),
    aspect.ratio=1,
    strip.text.y=element_blank(),
    legend.position="top"
  )
  gp = gp + xlab("Low concentration log2(fold change)")
  gp = gp + ylab("High concentration log2(fold change)")
  if (!x_title) {
    gp = gp + theme(axis.title.x=element_blank())
  }
  if (!y_title) {
    gp = gp + ylab("")
  }
  gp
}

library(ggpubr)

outfile = file.path(outdir, "fc_correlation.png")

met_split = split(metabolites_2, ceiling(seq_along(metabolites_2)/9))

png(outfile, width=1800, height=2200, res=250)
ggarrange(
  make_plot(met_split[[1]], F, F),
  make_plot(met_split[[2]], F, T),
  make_plot(met_split[[3]], T, F),
  common.legend=T,
  nrow = 3,
  ncol = 1,
  heights=unlist(lapply(met_split, length)) * c(1,0.79,1.087)
)
dev.off()
