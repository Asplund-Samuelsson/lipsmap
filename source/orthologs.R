options(width=110)
library(tidyverse)

# Load data
lipsmap = read_tsv(
  "data/annotated_comparison_results.tab.gz",
  col_types = cols(Protein.names = col_character())
)

orthologs = read_tsv(
  "data/uniprot_eggNOG.tab",
  col_names = c("UniProt_entry", "Ortholog")
)

org_uniprot = read_tsv(
  "data/organism_uniprot.tab",
  col_names = c("Organism", "UniProt_entry")
)

# Select orthologous groups with more than one member
orthologs = orthologs %>%
  # Count proteins per ortholog group
  group_by(Ortholog) %>%
  summarise(Count = length(UniProt_entry)) %>%
  # Select only ortholog groups with more than one protein member
  filter(Count > 1) %>%
  select(Ortholog) %>%
  # Use selected orthologs to filter the full ortholog table
  inner_join(orthologs)

# Determine ortholog interactions
ortholog_interactions = lipsmap %>%
  group_by(Organism, Metabolite, Conc, UniProt_entry) %>%
  summarise(Interaction = "sign" %in% Sign) %>%
  inner_join(orthologs) %>%
  mutate(Conc = ifelse(Conc == 1, "Low", "High"))

# Simplify and filter
interaction_comparison = ortholog_interactions %>%
  # Select only High Conc
  filter(Conc == "High") %>%
  # Simplify to existence of interaction within ortholog
  group_by(Organism, Metabolite, Ortholog) %>%
  summarise(Interaction = TRUE %in% Interaction) %>%
  # Keep only data with at least one Interaction and more than one example
  group_by(Ortholog) %>%
  filter(length(Interaction) > 1, sum(Interaction) > 0) %>%
  # Sort Orthologs and metabolites by number of Interaction
  ungroup() %>%
  mutate(
    Metabolite = factor(
      Metabolite,
      levels = (
        group_by(., Metabolite) %>%
        summarise(Count = sum(Interaction)) %>%
        arrange(-Count, Metabolite) %>%
        pull(Metabolite)
      )
    ),
    Ortholog = factor(
      Ortholog,
      levels = (
        group_by(., Ortholog) %>%
        summarise(Count = sum(Interaction)) %>%
        arrange(-Count, Ortholog) %>%
        pull(Ortholog)
      )
    )
  )

# Plot tiles
gp = ggplot(
  interaction_comparison,
  aes(x=Ortholog, y=Organism, fill=Interaction)
)
gp = gp + geom_tile()
gp = gp + facet_grid(Metabolite~.)
gp = gp + theme_bw()
gp = gp + scale_fill_manual(values=c("#9970ab", "#a6dba0"))
gp = gp + theme(
  axis.ticks = element_line(colour="black"),
  axis.ticks.x = element_blank(),
  axis.text.x = element_blank(),
  axis.text = element_text(colour="black"),
  strip.background = element_blank(),
  panel.grid = element_blank(),
  strip.text.y = element_text(angle=0, hjust=0)
)

ggsave("results/orthologs_interaction_comparison.png", gp, w=16, h=10)

# Calculate Jaccard similarity
jaccard = function(x, y){
  a = sum(x == 1 & y == 1)
  b = sum(x != y)
  return(a / (a + b))
}

# Jaccard similarity disregards cases where both vectors are zero
# ...therefore it can be used directly on the comparison data already acquired
interaction_comparison_wide = interaction_comparison %>%
  group_by(Organism, Metabolite) %>%
  spread(Ortholog, Interaction) %>%
  ungroup()

jaccard_pairs = as_tibble(t(combn(1:nrow(interaction_comparison_wide), 2)))

library(foreach)
library(doMC)
registerDoMC(16)

ortholog_jaccard = bind_rows(foreach(i=1:nrow(jaccard_pairs)) %dopar% {
  # Select data for pair to compare
  a = interaction_comparison_wide[jaccard_pairs[i,]$V1,]
  b = interaction_comparison_wide[jaccard_pairs[i,]$V2,]
  # Extract binary vectors
  a_bin = a %>% select(-Organism, -Metabolite) %>% as.numeric()
  b_bin = b %>% select(-Organism, -Metabolite) %>% as.numeric()
  # Subset to elements that are not NA; orthologs must exist
  not_missing = !(is.na(a_bin) | is.na(b_bin))
  j = jaccard(a_bin[not_missing], b_bin[not_missing])
  # Create table with results
  tibble(
    Organism_A = a$Organism, Metabolite_A = a$Metabolite,
    Interactions_A = sum(a_bin[not_missing]),
    Organism_B = b$Organism, Metabolite_B = b$Metabolite,
    Interactions_B = sum(b_bin[not_missing]),
    Jaccard = j, Orthologs = sum(not_missing)
  )
})

# Perform Principal Coordinates Analysis
library(vegan)

ortholog_jaccard_dist = interaction_comparison_wide %>%
  select(-Organism, -Metabolite) %>%
  as.matrix() %>%
  vegdist("jaccard", na.rm=T)

ortholog_jaccard_dist_m = as.matrix(ortholog_jaccard_dist)
ortholog_jaccard_dist_m[!is.finite(ortholog_jaccard_dist_m)] = 1
ortholog_jaccard_dist = as.dist(ortholog_jaccard_dist_m)

library(ape)

ortholog_jaccard_pcoa = pcoa(ortholog_jaccard_dist)

# Create plotting table
ortholog_jaccard_pcoa_plot = ortholog_jaccard_pcoa$vectors[,1:2] %>%
  as_tibble() %>%
  mutate(
    Organism = interaction_comparison_wide$Organism,
    Metabolite = interaction_comparison_wide$Metabolite
  ) %>%
  rename(PCo1 = Axis.1, PCo2 = Axis.2) %>%
  inner_join(
    interaction_comparison %>%
      group_by(Organism, Metabolite) %>%
      summarise(Interactions = sum(Interaction))
  )

# Calculate fraction of variance per PC
library(scales)
ortholog_jaccard_pcoa_var = percent(ortholog_jaccard_pcoa$values$Relative_eig)

# Plot it
library(ggrepel)

gp = ggplot(
  ortholog_jaccard_pcoa_plot,
  aes(x=PCo1, y=PCo2, colour=Organism, label=Metabolite)
)
gp = gp + geom_point(alpha=0.8, mapping=aes(size=Interactions))
gp = gp + scale_colour_manual(values=c("#9970ab","#a6dba0","#1b7837"), guide=F)
gp = gp + labs(
            x=paste("PCo1 (", ortholog_jaccard_pcoa_var[1], ")", sep=""),
            y=paste("PCo2 (", ortholog_jaccard_pcoa_var[2], ")", sep="")
          )
gp = gp + theme_bw()
gp = gp + facet_wrap(~Metabolite, ncol=5)
gp = gp + theme(
  axis.ticks = element_line(colour="black"),
  axis.text = element_text(colour="black"),
  strip.background = element_blank(),
  aspect.ratio = 1
)
gp1 = gp

gp = ggplot(
  ortholog_jaccard_pcoa_plot,
  aes(x=PCo1, y=PCo2, colour=Organism, label=Metabolite)
)
gp = gp + geom_point(alpha=0.8, mapping=aes(size=Interactions))
gp = gp + geom_text_repel(force=3, size=3,alpha=0.9)
gp = gp + scale_colour_manual(values=c("#9970ab","#a6dba0","#1b7837"), guide=F)
gp = gp + labs(
            x=paste("PCo1 (", ortholog_jaccard_pcoa_var[1], ")", sep=""),
            y=paste("PCo2 (", ortholog_jaccard_pcoa_var[2], ")", sep="")
          )
gp = gp + theme_bw()
gp = gp + facet_wrap(~Organism, ncol=1)
gp = gp + theme(
  axis.ticks = element_line(colour="black"),
  axis.text = element_text(colour="black"),
  strip.background = element_blank(),
  aspect.ratio = 1
)
gp2 = gp

library(ggpubr)
outfile = "results/orthologs_interaction_pcoa.pdf"
pdf(outfile, width=11, height=8.5, onefile=FALSE)
ggarrange(
  gp2,
  gp1,
  nrow = 1,
  ncol = 2,
  labels = c("A", "B"),
  widths = c(1, 2),
  common.legend = T,
  label.x = 0, label.y = 1
)
garbage = dev.off()

# Save Jaccard index table
write_tsv(
  arrange(ortholog_jaccard, -Jaccard),
  "results/orthologs_interaction_jaccard.tab"
)

# Determine orthologs that are always detected
omnipresent_orthologs = interaction_comparison_wide %>%
  select(-Organism, -Metabolite) %>%
  as.matrix() %>%
  colSums()

omnipresent_orthologs = omnipresent_orthologs[!is.na(omnipresent_orthologs)]

# Run PCA on raw orthologs interaction data
ortholog_pca_data = interaction_comparison_wide #%>%
#  select(Organism, Metabolite, all_of(names(omnipresent_orthologs)))

ortholog_pca_data_m = ortholog_pca_data %>%
  select(-Organism, -Metabolite) %>%
  data.matrix() %>%
  replace_na(0)

ortholog_pca = prcomp(ortholog_pca_data_m)

# Create plotting dataframes
ortholog_pca_plot = as_tibble(ortholog_pca$x) %>% select(PC1, PC2)
ortholog_pca_plot = bind_cols(
  ortholog_pca_data %>% select(Organism, Metabolite),
  ortholog_pca_plot,
  tibble(Interactions = rowSums(ortholog_pca_data_m))
)

# Calculate fraction of variance per PC
library(scales)
ortholog_pca_var = percent(ortholog_pca$sdev^2 / sum(ortholog_pca$sdev^2))

# Plot it
gp = ggplot(
  ortholog_pca_plot,
  aes(x=PC1, y=PC2, colour=Organism, label=Metabolite)
)
gp = gp + geom_point(alpha=0.8, mapping=aes(size=Interactions))
gp = gp + scale_colour_manual(values=c("#9970ab","#a6dba0","#1b7837"), guide=F)
gp = gp + scale_size_continuous(breaks=c(0,20,100,200))
gp = gp + labs(
            x=paste("PC1 (", ortholog_pca_var[1], ")", sep=""),
            y=paste("PC2 (", ortholog_pca_var[2], ")", sep="")
          )
gp = gp + theme_bw()
gp = gp + facet_wrap(~Metabolite, ncol=5)
gp = gp + theme(
  axis.ticks = element_line(colour="black"),
  axis.text = element_text(colour="black"),
  strip.background = element_blank(),
  aspect.ratio = 1
)
gp1 = gp

gp = ggplot(
  ortholog_pca_plot,
  aes(x=PC1, y=PC2, colour=Organism, label=Metabolite)
)
gp = gp + geom_point(alpha=0.8, mapping=aes(size=Interactions))
gp = gp + geom_text_repel(force=3, size=3,alpha=0.9)
gp = gp + scale_size_continuous(breaks=c(0,20,100,200))
gp = gp + scale_colour_manual(values=c("#9970ab","#a6dba0","#1b7837"), guide=F)
gp = gp + labs(
            x=paste("PC1 (", ortholog_pca_var[1], ")", sep=""),
            y=paste("PC2 (", ortholog_pca_var[2], ")", sep="")
          )
gp = gp + theme_bw()
gp = gp + facet_wrap(~Organism, ncol=1)
gp = gp + theme(
  axis.ticks = element_line(colour="black"),
  axis.text = element_text(colour="black"),
  strip.background = element_blank(),
  aspect.ratio = 1
)
gp2 = gp

library(ggpubr)
outfile = "results/orthologs_interaction_pca.pdf"
pdf(outfile, width=11, height=8.5, onefile=FALSE)
ggarrange(
  gp2,
  gp1,
  nrow = 1,
  ncol = 2,
  labels = c("A", "B"),
  widths = c(1, 2),
  common.legend = T,
  label.x = 0, label.y = 1
)
garbage = dev.off()
