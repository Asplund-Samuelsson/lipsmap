options(width=110)
library(ggtree)
library(vegan)
library(phytools)
library(ape)
library(scales)
library(ggrepel)
library(tidyverse)
library(ggpubr)

# Read input file and output directory
args = commandArgs(trailingOnly=TRUE)
infile = args[1]
outdir = args[2]

# Load data
lipsmap = read_tsv(infile) %>%
  # Change significance
  mutate(Sign = ifelse(adj.pvalue < 0.01, "sign", "unsign"))

# Define organisms and colors
organisms = c("Hydrogenophaga", "Cupriavidus", "Synechococcus", "Synechocystis")
organcols = c("#762a83", "#9970ab","#5aae61","#1b7837")
names(organcols) = organisms

# Determine interactions
interactions = lipsmap %>%
  # Determine metabolite-interacting proteins; at least one significant peptide
  group_by(Organism, Metabolite, Conc, UniProt_entry) %>%
  summarise(Interaction = "sign" %in% Sign) %>%
  # Clarify concentration description
  mutate(Conc = ifelse(Conc == 1, "Low", "High"))

# Determine what experiments were done
experiments = interactions %>%
  ungroup() %>%
  select(-Interaction, -UniProt_entry) %>%
  distinct()

# Determine what proteins were detected at least once
proteins = interactions %>%
  ungroup() %>%
  select(Organism, UniProt_entry) %>%
  distinct()

# Determine possible detected proteins in all experiments
possible = inner_join(experiments, proteins)

# Filter out proteins not detected in all experiments
filtered_interactions = interactions %>%
  anti_join(
    possible %>%
      left_join(interactions) %>%
      filter(is.na(Interaction)) %>%
      select(UniProt_entry) %>%
      distinct()
  )

# Function to get PCA results
do_pca = function (tb) {
  # Spread interactions into columns by protein
  tb_wide = tb %>%
    ungroup() %>%
    select(Metabolite, UniProt_entry, Interaction) %>%
    spread(UniProt_entry, Interaction)

  # Make matrix
  tb_matrix = tb_wide %>% select(-Metabolite) %>% as.matrix() * 1

  # Do PCA
  tb_pca = prcomp(tb_matrix)

  # Calculate variance
  tb_pca_var = percent(tb_pca$sdev^2 / sum(tb_pca$sdev^2))

  # Return data
  tb_wide %>%
    select(Metabolite) %>%
    bind_cols(as_tibble(tb_pca$x[,1:2])) %>%
    gather(PC, Value, -Metabolite) %>%
    inner_join(tibble(PC = c("PC1", "PC2"), Variance = tb_pca_var[1:2])) %>%
    inner_join(tb %>% select(Organism, Metabolite, Conc) %>% distinct())
}

# Do all PCA
pca_results = filtered_interactions %>%
  # Do PCA on each organism and concentration
  group_by(Organism, Conc) %>%
  group_split() %>%
  lapply(do_pca) %>%
  bind_rows() %>%
  # Normalize PCs between -1 and 1
  group_by(Organism, Conc, PC) %>%
  mutate(
    Norm = -1 + (Value - min(Value))*(1-(-1))/(max(Value) - min(Value))
  ) %>%
  # Order organisms
  mutate(
    Organism = factor(Organism, levels=organisms),
    Conc = paste(Conc, "concentration", sep=" ")
  )

# Create variance table
pca_var = pca_results %>%
  ungroup() %>%
  select(PC, Variance, Organism, Conc) %>%
  distinct() %>%
  mutate(
    PC1 = ifelse(PC == "PC1", 0.8, -0.8),
    PC2 = ifelse(PC == "PC1", -1.1, 1.1),
    Variance = paste(PC, ": ", Variance, sep="")
  )

# Select PCA plotting data
pca_plot = pca_results %>%
  ungroup() %>%
  select(Organism, Conc, Metabolite, PC, Norm) %>%
  spread(PC, Norm) %>%
  # Keep only interesting metabolites
  mutate(
    Metabolite = ifelse(
      Metabolite %in% c("NADPH", "NADP", "ADP", "AMP", "G6P", "FBP"),
      Metabolite, NA
    ),
    Show = !is.na(Metabolite)
  )

# Plot it
gp = ggplot(
  pca_plot, aes(color=Organism, label=Metabolite, x=PC1, y=PC2, alpha=Show)
)
gp = gp + geom_point()
gp = gp + geom_text_repel(
  max.overlaps = 1000
)
gp = gp + geom_label(
  data=pca_var, mapping=aes(label=Variance), color="#666666", alpha=0.3, size=3
)
gp = gp + scale_color_manual(values=organcols, guide="none")
gp = gp + scale_alpha_manual(values=c(0.3, 1), guide="none")
gp = gp + facet_grid(Conc ~ Organism)
gp = gp + xlim(-1.2, 1.2) + ylim(-1.2, 1.2)
gp = gp + theme_bw()
gp = gp + theme(
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  panel.grid = element_blank(),
  aspect.ratio=1,
  strip.background=element_blank()
)

ggsave(file.path(outdir, "metabolites_pca.pdf"), gp, w=10, h=5)

# Function to get hierarchical clustering results
do_clust = function (tb) {

  # Spread interactions into columns by protein
  tb_wide = tb %>%
    ungroup() %>%
    select(Metabolite, UniProt_entry, Interaction) %>%
    spread(UniProt_entry, Interaction)

  # Make matrix
  tb_matrix = tb_wide %>% select(-Metabolite) %>% as.matrix() * 1

  # Add rownames (metabolites)
  rownames(tb_matrix) = pull(tb_wide, Metabolite)

  # Do hierarchical clustering
  tb_dist = vegdist(tb_matrix, "jaccard")

  # Replace missing values with maximum distance (1)
  tb_dist_m = as.matrix(tb_dist)
  tb_dist_m[!is.finite(tb_dist_m)] = 1
  tb_dist = as.dist(tb_dist_m)

  # Cluster based on distance
  tb_clust = hclust(tb_dist, method = "ward.D2")

  # Count interactions
  tb_count = tb %>%
    group_by(Metabolite) %>%
    summarise(Interactions = sum(Interaction)) %>%
    rename(label = Metabolite)

  # Create tree from clustering
  tb_tree = as.phylo(tb_clust)

  # Plot tree
  gp = ggtree(
    tb_tree, aes(colour=Organism, fill=Organism), layout="rectangular"
  )

  # Add Organism to plot data
  gp$data = gp$data %>% mutate(Organism = unique(tb$Organism))

  # Configure tree
  gp = gp + geom_tiplab(aes(label=label), size=2.5)
  gp = gp + scale_color_manual(values=organcols, guide="none")
  gp = gp + scale_fill_manual(values=organcols, guide="none")
  gp = gp + scale_x_continuous(
    expand = expansion(add = c(0.05, 0.5)), n.breaks=4
  )

  # Add bar plot with interaction counts
  gp = gp + geom_facet(
    panel = 'Interactions',
    data = tb_count %>% select(label, Interactions),
    geom = geom_col,
    mapping = aes(x=Interactions),
    orientation = 'y'
  )

  gp = gp + theme_tree2()
  gp = gp + theme(
    strip.background = element_blank(),
    strip.text = element_blank(),
    plot.margin = unit(c(0.3,0.3,0.3,0.3), "cm"),
    axis.text.x = element_text(size=8)
  )

  # Return finished tree plot
  return(gp)
}

# Do all hierarchical clustering
clust_groups = filtered_interactions %>%
  # Do PCA on each organism and concentration
  mutate(
    Organism = factor(Organism, levels=organisms),
    Conc = factor(
      paste(Conc, " concentration", sep=""),
      levels = c("High concentration", "Low concentration")
    )
  ) %>%
  group_by(Conc, Organism) %>%
  group_split()

# Pre-allocate list
clust_results = vector("list", length(clust_groups))

for (x in 1:length(clust_groups)) {
  clust_results[[x]] = do_clust(clust_groups[[x]])
}

# Make label column and header
make_title = function(title_text, axis) {
  gp = ggplot(data.frame())
  gp = gp + geom_blank()
  gp = gp + theme_void()
  if (axis == "x"){
    gp = gp + xlab(title_text)
    gp = gp + scale_x_continuous(position = "top")
  } else {
    gp = gp + ylab(title_text)
    gp = gp + scale_y_continuous(position = "right")
    gp = gp + theme(axis.title=element_text(angle=270))
  }
  gp = gp + theme(axis.title=element_text(color="black"))
  return(gp)
}

# Create header with column titles
plot_header = c(
  lapply(
    unique(unlist(lapply(
      clust_groups, function(x){unique(pull(x, Organism))}
    ))),
    make_title,
    axis="x"
  ),
  list(NULL)
)

# Create column with row titles
plot_rows = lapply(
  unique(unlist(lapply(
    clust_groups, function(x){unique(pull(x, Conc))}
  ))),
  make_title,
  axis="y"
)

# Make list of all plots in correct order
plot_list = c(
  plot_header,
  unlist(lapply(
    1:length(plot_rows),
    function(x){c(
      # Split actual plots into rows for each concentration
      split(
        clust_results,
        cut(
          seq_along(clust_results),
          length(plot_rows),
          labels = FALSE
        )
      )[[x]],
      list(plot_rows[[x]])
    )}
  ), recursive = F)
)

# Plot it
outfile = file.path(outdir, "metabolites_jaccard_ward_D2.pdf")
pdf(outfile, width=10, height=7, onefile=FALSE)
ggarrange(
  plotlist=plot_list,
  common.legend=T,
  nrow = length(plot_rows) + 1,
  ncol = length(plot_header),
  widths = c(rep(1, length(plot_header)-1), 0.1),
  heights = c(0.05, rep(1, length(plot_rows)))
)
dev.off()
