options(width=110)
library(tidyverse)
library(scales)

# Load data
lipsmap = bind_rows(
  read_csv("data/accoa_20210127.csv") %>% mutate(Date = "Old"),
  read_csv("data/accoa_20210527.csv") %>% mutate(Date = "New")
) %>%
  select(Peptide_ID, Conc, adj.pvalue, log2FC, Date, UniProt_entry) %>%
  mutate(Conc = ifelse(Conc == 1, "Low", "High"))

# Correlate q value between Old and New
q_comparison = lipsmap %>%
  group_by(Conc, Date, UniProt_entry) %>%
  summarise(q = min(adj.pvalue)) %>%
  spread(Date, q)

# Low correlation
r_low = cor(
  filter(q_comparison, Conc == "Low")$New,
  filter(q_comparison, Conc == "Low")$Old,
  use="complete", method="spearman"
)

# High correlation
r_high = cor(
  filter(q_comparison, Conc == "High")$New,
  filter(q_comparison, Conc == "High")$Old,
  use="complete", method="spearman"
)

# Full correlation
r_full = cor(
  q_comparison$New,
  q_comparison$Old,
  use="complete", method="spearman"
)

# Determine limits
q_lims = c(
  10**floor(log10(min(c(q_comparison$Old, q_comparison$New), na.rm=T))),
  10**ceiling(log10(max(c(q_comparison$Old, q_comparison$New), na.rm=T)))
)

# Plot q value comparison
gp = ggplot(
  mutate(q_comparison, Conc = paste(Conc, " concentration", sep="")),
  aes(x=Old, y=New)
)
gp = gp + geom_point(alpha=0.5, stroke=0, size=1, colour="#1b7837")
gp = gp + geom_text(
  data=tibble(
    Label = c(
      paste("R = ", round(r_low, 3), sep=""),
      paste("R = ", round(r_high, 3), sep="")
    ),
    Conc = c("Low concentration", "High concentration"),
    Old = c(1e-6, 1e-6),
    New = c(5e-7, 5e-7)
  ),
  mapping=aes(label=Label)
)
gp = gp + theme_bw()
gp = gp + facet_grid(.~Conc)
gp = gp + scale_x_log10(limits=q_lims)
gp = gp + scale_y_log10(limits=q_lims)
gp = gp + theme(
  axis.text = element_text(color="black"),
  axis.ticks = element_line(color="black"),
  strip.background = element_blank(),
  aspect.ratio = 1
)
gp = gp + xlab("min(q) on 2021-01-27")
gp = gp + ylab("min(q) on 2021-05-27")
gp = gp + ggtitle(
  paste(
    "Correlation of min(q) between dates (overall R = ",
    round(r_full, 3), ")", sep=""
  )
)

ggsave("results/accoa_min_q_cor.pdf", gp, h=3.5, w=7)

# Investigate interactions
interaction = q_comparison %>%
  # Add some different cutoff alternatives
  uncount(5) %>%
  group_by(Conc, UniProt_entry) %>%
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
  group_by(Conc, Cutoff, Interaction) %>%
  summarise(Proteins = length(UniProt_entry)) %>%
  # Order Interaction class
  mutate(
    Interaction = factor(Interaction, levels=c("Yes", "Split", "No"))
  ) %>%
  # Add missing values
  complete(Interaction, fill=list(Proteins = 0)) %>%
  # Calculate percent of each Interaction type
  group_by(Conc, Cutoff) %>%
  mutate(
    Percent = percent(Proteins / sum(Proteins)),
    Percent = ifelse(Proteins / sum(Proteins) < 0.1, "", Percent),
    Conc = paste(Conc, "concentration", sep=" ")
  )

# Plot it
gp = ggplot(
  interaction,
  aes(x=Cutoff, y=Proteins, fill=Interaction, label=Percent)
)
gp = gp + geom_col(color="black", position=position_dodge(width=0.7), width=0.7)
gp = gp + geom_text(
  position=position_dodge(width=0.7), angle=90, hjust=1.1, vjust=0.5, size=2.5
)
gp = gp + theme_bw()
gp = gp + facet_grid(.~Conc)
gp = gp + scale_x_log10()
gp = gp + scale_fill_brewer(palette="PuOr")
gp = gp + theme(
  axis.text = element_text(color="black"),
  axis.ticks = element_line(color="black"),
  strip.background = element_blank()
)
gp = gp + xlab("q value cutoff")

ggsave("results/accoa_date_agreement.pdf", gp, h=3, w=7)
