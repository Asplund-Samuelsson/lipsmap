options(width=110)
library(tidyverse)

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
interactions = lipsmap %>%
  full_join(tibble(Threshold = c(0.01, 0.001, 0.0001))) %>%
  group_by(Date, UniProt_entry, Conc) %>%
  # Interacting proteins have
  summarise(Interaction = sum(adj.pvalue < 0.01) > 0)
