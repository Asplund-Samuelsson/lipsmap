options(width=110)
library(tidyverse)
library(scales)

# Load data
lipsmap = bind_rows(
  read_csv("data/accoa_20210127.csv") %>%
    mutate(Metabolite = "AcCoA", Date = "Old"),
  read_csv("data/accoa_20210527.csv") %>%
    mutate(Metabolite = "AcCoA", Date = "New"),
  read_csv("data/atp_20200612.csv") %>%
    mutate(Metabolite = "ATP", Date = "Old"),
  read_csv("data/atp_20200819.csv") %>%
    mutate(Metabolite = "ATP", Date = "New")
) %>%
  select(
    Peptide_ID, Conc, adj.pvalue, log2FC,
    Date, Metabolite, UniProt_entry
  ) %>%
  mutate(
    Conc = case_when(
      Conc == 2 ~ "High",
      Conc == 1 & Metabolite == "ATP" ~ "High",
      T ~ "Low"
    )
  )

# Calculate top 25% by number of peptides per protein
top_proteins = lipsmap %>%
  group_by(Metabolite, Conc, Date, UniProt_entry) %>%
  summarise(Peptides = length(Peptide_ID)) %>%
  top_frac(0.25, Peptides)

# Classify proteins as being in the top 25% or not
lipsmap = bind_rows(
  lipsmap %>%
    mutate(Peptides = "All proteins"),
  lipsmap %>%
    inner_join(top_proteins) %>%
    mutate(Peptides = "Top 25% peptide count")
)

# Correlate q value between Old and New
q_comparison = lipsmap %>%
  group_by(Metabolite, Conc, Date, UniProt_entry, Peptides) %>%
  summarise(q = min(adj.pvalue)) %>%
  spread(Date, q) %>%
  mutate(Condition = paste(Metabolite, " (", Conc, " concentration)", sep=""))

# Determine limits
q_lims = c(
  10**floor(log10(min(c(q_comparison$Old, q_comparison$New), na.rm=T))),
  10**ceiling(log10(max(c(q_comparison$Old, q_comparison$New), na.rm=T)))
)

# Correlate data
q_correlation = q_comparison %>%
  group_by(Metabolite, Conc, Peptides) %>%
  summarise(Correlation = cor(New, Old, use="complete", method="spearman")) %>%
  left_join(
    q_comparison %>%
      ungroup() %>%
      select(Metabolite, Conc, Condition) %>%
      distinct()
  ) %>%
  ungroup() %>%
  mutate(
    Old = rep(q_lims[1]*10, length(Correlation)),
    New = rep(q_lims[1]*5, length(Correlation)),
    Label = paste("R = ", round(Correlation, 3), sep=""),
  )

# Plot q value comparison
gp = ggplot(q_comparison, aes(x=Old, y=New))
gp = gp + geom_point(mapping=aes(colour=Condition), alpha=0.5, stroke=0, size=1)
gp = gp + scale_color_manual(values=c("#1b7837", "#1b7837", "#9970ab"), guide=F)
gp = gp + geom_text(data=q_correlation, mapping=aes(label=Label))
gp = gp + theme_bw()
gp = gp + facet_grid(Peptides~Condition)
gp = gp + scale_x_log10(limits=q_lims)
gp = gp + scale_y_log10(limits=q_lims)
gp = gp + theme(
  axis.text = element_text(color="black"),
  axis.ticks = element_line(color="black"),
  strip.background = element_blank(),
  aspect.ratio = 1
)
gp = gp + xlab("min(q) in early experiment")
gp = gp + ylab("min(q) in later experiment")

ggsave("results/dates_min_q_cor.pdf", gp, h=5.5, w=8)

# Investigate interactions
interaction = q_comparison %>%
  # Add some different cutoff alternatives
  uncount(5) %>%
  group_by(Condition, UniProt_entry, Peptides) %>%
  mutate(Cutoff = 10**c(-5:-1)) %>%
  ungroup() %>%
  # Determine what the two dates are agreeing upon
  mutate(
    Interaction = case_when(
      is.na(New) | is.na(Old) ~ "Single",
      New < Cutoff & Old < Cutoff ~ "Yes",
      New >= Cutoff & Old >= Cutoff ~ "No",
      T ~ "Split"
    )
  ) %>%
  # Remove single proteins (detected in only one experiment)
  filter(Interaction != "Single") %>%
  # Summarise Interaction decision
  group_by(Condition, Cutoff, Interaction, Peptides) %>%
  summarise(Proteins = length(UniProt_entry)) %>%
  # Order Interaction class
  mutate(
    Interaction = factor(Interaction, levels=c("Yes", "Split", "No"))
  ) %>%
  # Add missing values
  ungroup() %>%
  complete(
    expand(., Peptides, Interaction, Condition, Cutoff),
    fill=list(Proteins = 0)
  ) %>%
  # Calculate percent of each Interaction type
  group_by(Condition, Cutoff, Peptides) %>%
  mutate(
    Percent = percent(Proteins / sum(Proteins)),
    Percent = ifelse(Proteins / sum(Proteins) < 0.1, "", Percent)
  )

# Determine frequency of interaction at different Conc and Cutoff
exp_int = q_comparison %>%
  # Add the cutoffs
  uncount(5) %>%
  group_by(Condition, UniProt_entry, Peptides) %>%
  mutate(Cutoff = 10**c(-5:-1)) %>%
  # Classify interaction
  mutate(New = New < Cutoff, Old = Old < Cutoff) %>%
  # Calculate fraction interaction
  gather(
    Date, Interaction,
    -Metabolite, -Conc, -Condition, -UniProt_entry, -Cutoff, -Peptides
  ) %>%
  filter(!is.na(Interaction)) %>%
  group_by(Condition, Date, Cutoff, Peptides) %>%
  summarise(Interaction = sum(Interaction)/length(Interaction)) %>%
  # Spread dates into columns
  spread(Date, Interaction) %>%
  # Add protein counts (only proteins that were detected on both dates)
  left_join(
    interaction %>%
      group_by(Condition, Cutoff, Peptides) %>%
      summarise(Proteins = sum(Proteins))
  ) %>%
  # Calculate expected agreements if randomly sampled
  mutate(
    Yes = New * Old * Proteins,
    No = (1-New) * (1-Old) * Proteins,
    Split = Proteins - Yes - No
  ) %>%
  # Gather into long format
  select(-Proteins, -New, -Old) %>%
  gather(Interaction, Proteins, -Condition, -Cutoff, -Peptides) %>%
  # Order Interaction
  mutate(Interaction = factor(Interaction, levels=c("Yes", "Split", "No")))

# Plot it
gp = ggplot(
  interaction,
  aes(x=Cutoff, y=Proteins, fill=Interaction)
)
gp = gp + geom_col(color="black", position=position_dodge(width=0.7), width=0.7)
gp = gp + geom_point(
  data=exp_int, shape=21, alpha=1, position=position_dodge(width=0.7),
  fill="white", colour="black", size=1, mapping=aes(group=Interaction)
)
gp = gp + geom_text(
  position=position_dodge(width=0.7), angle=90, hjust=1.1, vjust=0.5, size=2.5,
  mapping=aes(label=Percent)
)
gp = gp + theme_bw()
gp = gp + facet_wrap(Peptides~Condition, scales="free_y")
gp = gp + scale_x_log10()
gp = gp + scale_fill_brewer(palette="PuOr")
gp = gp + theme(
  axis.text = element_text(color="black"),
  axis.ticks = element_line(color="black"),
  strip.background = element_blank(),
  aspect.ratio = 1,
  legend.position = "bottom"
)
gp = gp + xlab("q value cutoff")
gp = gp + ylab("Proteins (circle = random)")

ggsave("results/dates_agreement.pdf", gp, h=8, w=10.5)
