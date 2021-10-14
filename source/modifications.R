options(width=110)
library(tidyverse)

# Load data
lipsmap = read_tsv("data/annotated_comparison_results.tab.gz") %>%
  # Change significance
  mutate(Sign = ifelse(adj.pvalue < 0.01, "sign", "unsign"))

uniprot_locus = read_tsv(
    "data/uniprot_locus_complete.tab", col_names=c("UniProt_entry", "Locus")
)

acetylation = scan("data/synechocystis_acetylation.txt", character())
propionylation = scan("data/synechocystis_propionylation.txt", character())

interaction = lipsmap %>%
# Consider interactions only in Synechocystis, with AcCoA and GAP
  filter(Organism == "Synechocystis", Metabolite %in% c("GAP", "AcCoA")) %>%
  # Determine significant proteins
  group_by(Metabolite, Conc, UniProt_entry) %>%
  summarise(Interaction = "sign" %in% Sign) %>%
  # Add locus
  inner_join(uniprot_locus) %>%
  mutate(
    # Clarify concentration
    Conc = ifelse(Conc == 1, "Low", "High"),
    # Determine acetylation status
    Acetylation = Locus %in% acetylation,
    # Determine propionylation status
    Propionylation = Locus %in% propionylation
  ) %>%
  # Gather the modifications
  select(-UniProt_entry) %>%
  gather(Modification, Modified, -Metabolite, -Conc, -Locus, -Interaction)

fisher_results = interaction %>%
  # Count groups for Fisher test
  group_by(Metabolite, Conc, Modification) %>%
  summarise(
    IM = sum(Interaction & Modified),
    iM = sum(!Interaction & Modified),
    Im = sum(Interaction & !Modified),
    im = sum(!Interaction & !Modified)
  ) %>%
  # Do Fisher test
  group_by(Metabolite, Conc, Modification) %>%
  mutate(
    p = fisher.test(matrix(c(IM, iM, Im, im), ncol=2, byrow=T))$p,
    Significant = p < 0.01
  )

# Prepare plotting data
plot_data = fisher_results %>%
  select(-p) %>%
  gather(Class, Proteins, -Metabolite, -Conc, -Modification, -Significant) %>%
  mutate(
    Interaction = grepl("I", Class),
    Modified = grepl("M", Class),
    Modification = str_replace(Modification, "ion$", "ed"),
    Modification_type = Modification,
    Modification = factor(
      ifelse(
        Modified, Modification,
        paste("Not ", str_to_lower(Modification), sep="")
      ),
      levels=c(
        "Acetylated", "Not acetylated", "Propionylated", "Not propionylated"
      )
    ),
    Conc = paste(Conc, " concentration", sep="")
  )

# Create percent label
library(scales)
percent_labels = plot_data %>%
  select(
    Metabolite, Conc, Modification, Modification_type,
    Proteins, Interaction, Modified
  ) %>%
  mutate(Interaction = ifelse(Interaction, "Yes", "No")) %>%
  spread(Interaction, Proteins) %>%
  mutate(
    Label = percent(round(Yes / (Yes + No), 3)),
    Proteins = Yes,
    Significant = T
  )

# Plot it
gp = ggplot(
  plot_data,
  aes(
    x=Modification, y=Proteins,
    fill=Interaction, color=Interaction,
    linetype=Significant, alpha=Modified
  )
)
gp = gp + geom_col(position="stack", size=0.5, width=0.9)
gp = gp + scale_fill_manual(values=c("#e08214", "#8073ac"))
gp = gp + scale_color_manual(values=c("#7f3b08", "#2d004b"))
gp = gp + geom_label(
  data=percent_labels, mapping=aes(label=Label), color="black", alpha=0.5,
  fill="white"
)
gp = gp + geom_rect(
    xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, colour="black",
    size=0.4, alpha=0.6
  )
gp = gp + scale_linetype_manual(values=c(2,1), guide=F)
gp = gp + scale_alpha_manual(
  values=c(0.6,1), guide = guide_legend(reverse = TRUE)
)
gp = gp + facet_grid(Metabolite ~ Conc + Modification_type, scales="free")
gp = gp + theme_bw()
gp = gp + theme(
  axis.text = element_text(color="black"),
  axis.ticks = element_line(color="black"),
  strip.background = element_blank(),
  legend.position = "top",
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.x = element_blank(),
  panel.border=element_blank()
)

ggsave("results/modifications.pdf", gp, w=7, h=7)
