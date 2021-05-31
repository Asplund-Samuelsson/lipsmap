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

eggnog_annotations = read_tsv("data/eggNOG_annotations.tab")

eggnog_category_descriptions = read_tsv(
  "data/ortholog_category_descriptions.tab"
)

# Determine eggNOG annotations
eggnog_annotations_unique = eggnog_annotations %>%
  # Select relevant data
  select(Ortholog, Version, Category) %>%
  # Split multiple categories between each character
  mutate(Category = str_split(Category, "")) %>%
  # Unnest categories per ortholog into rows
  unnest(cols="Category") %>%
  # Add missing version string to non-COG orthologs
  mutate(
    Ortholog = ifelse(
      str_starts(Ortholog, "COG"),
      Ortholog,
      ifelse(
        Version == 4.1,
        paste("ENOG41", Ortholog, sep=""),
        paste("ENOG50", Ortholog, sep="")
      )
    )
  ) %>%
  # Remove version and make rows unique
  select(-Version) %>%
  distinct()

# Determine ortholog interactions
ortholog_interactions = lipsmap %>%
  # Determine metabolite-interacting proteins; at least one significant peptide
  group_by(Organism, Metabolite, Conc, UniProt_entry) %>%
  summarise(Interaction = "sign" %in% Sign) %>%
  # Add ortholog information
  inner_join(orthologs) %>%
  # Clarify concentration description
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
  # Sum of matching interactions
  a = sum(x == 1 & y == 1)
  # Sum of mismatching interactions
  b = sum(x != y)
  # Return Jaccard similarity; matching interactions / total interactions
  return(a / (a + b))
}

# Jaccard similarity disregards cases where both vectors are zero
# ...therefore it can be used directly on the comparison data already acquired

# Spread interactions data into wide format to allow making matrices
interaction_comparison_wide = interaction_comparison %>%
  group_by(Organism, Metabolite) %>%
  spread(Ortholog, Interaction) %>%
  ungroup()

# Define non-redundant pairs to compare through Jaccard similarity
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
    # First organism and metabolite of pair; store number of interactions
    Organism_A = a$Organism, Metabolite_A = a$Metabolite,
    Interactions_A = sum(a_bin[not_missing]),
    # Second organism and metabolite of pair; store number of interactions
    Organism_B = b$Organism, Metabolite_B = b$Metabolite,
    Interactions_B = sum(b_bin[not_missing]),
    # Jaccard similarity and number of orthologs in comparison
    Jaccard = j, Orthologs = sum(not_missing)
  )
})

# Perform Principal Coordinates Analysis
library(vegan)

# Create Jaccard distance object
ortholog_jaccard_dist = interaction_comparison_wide %>%
  select(-Organism, -Metabolite) %>%
  as.matrix() %>%
  vegdist("jaccard", na.rm=T)

# Replace missing values with maximum distance (1)
ortholog_jaccard_dist_m = as.matrix(ortholog_jaccard_dist)
ortholog_jaccard_dist_m[!is.finite(ortholog_jaccard_dist_m)] = 1
ortholog_jaccard_dist = as.dist(ortholog_jaccard_dist_m)

library(ape)

ortholog_jaccard_pcoa = pcoa(ortholog_jaccard_dist)

# Create plotting table
ortholog_jaccard_pcoa_plot = ortholog_jaccard_pcoa$vectors[,1:2] %>%
  as_tibble() %>%
  # Add Organism and Metabolite
  mutate(
    Organism = interaction_comparison_wide$Organism,
    Metabolite = interaction_comparison_wide$Metabolite
  ) %>%
  # Give sensible names to axes
  rename(PCo1 = Axis.1, PCo2 = Axis.2) %>%
  # Add number of interactions per Organism and Metabolite
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

# Plot split by Metabolite
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

# Plot split by Organism
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

# Run PCA on raw orthologs interaction data

# Create matrix from interaction comparison data
interaction_comparison_wide_m = interaction_comparison_wide %>%
  select(-Organism, -Metabolite) %>%
  # Make binary
  data.matrix() %>%
  # Set missing values to zero
  replace_na(0)

ortholog_pca = prcomp(interaction_comparison_wide_m)

# Create plotting dataframes
ortholog_pca_plot = as_tibble(ortholog_pca$x) %>% select(PC1, PC2)
ortholog_pca_plot = bind_cols(
  interaction_comparison_wide %>% select(Organism, Metabolite),
  ortholog_pca_plot,
  tibble(Interactions = rowSums(interaction_comparison_wide_m))
)

# Calculate fraction of variance per PC
library(scales)
ortholog_pca_var = percent(ortholog_pca$sdev^2 / sum(ortholog_pca$sdev^2))

# Make plot split by Metabolite
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

# Make plot split by Organism
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

# Cluster orthologs
library(ggtree)
library(phytools)

# Simplify and filter
interaction_comparison_concs = ortholog_interactions %>%
  # Simplify to existence of interaction within ortholog
  group_by(Organism, Metabolite, Conc, Ortholog) %>%
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

# Make wide version
interaction_comparison_concs_wide = interaction_comparison_concs %>%
  group_by(Organism, Metabolite, Conc) %>%
  spread(Ortholog, Interaction) %>%
  ungroup()


# Make Jaccard distance object
ortholog_jaccard_dist_t = interaction_comparison_concs_wide %>%
  select(-Organism, -Metabolite, -Conc) %>%
  as.matrix() %>%
  # Transpose to cluster Orthologs
  t() %>%
  # Calculate Jaccard distance, removing positions with NA
  vegdist("jaccard", na.rm=T)

# Replace missing values with maximum distance (1)
ortholog_jaccard_dist_t_m = as.matrix(ortholog_jaccard_dist_t)
ortholog_jaccard_dist_t_m[!is.finite(ortholog_jaccard_dist_t_m)] = 1
ortholog_jaccard_dist_t = as.dist(ortholog_jaccard_dist_t_m)

# Create cluster object from distance object using ward.D2 method
ortholog_clustering = hclust(ortholog_jaccard_dist_t, method = "ward.D2")

# Divide into clusters
ortholog_clusters = cutree(ortholog_clustering, k=7)

# List clusters and Orthologs
ortholog_clusters = tibble(
  Ortholog = names(ortholog_clusters),
  Cluster = as.character(ortholog_clusters),
  label = Ortholog
)

# Make tree object from cluster object
ortholog_tree = as.phylo(ortholog_clustering)

# Find most recent common ancestor for each cluster, then colour by cluster
ortholog_mrca = ortholog_clusters %>%
  group_by(Cluster) %>%
  summarise(MRCA = findMRCA(ortholog_tree, label)) %>%
  mutate(node = lapply(MRCA, function(x){getDescendants(ortholog_tree, x)})) %>%
  unnest(cols="node")

# Plot clustering
gp = ggtree(
  as.phylo(ortholog_clustering),
  aes(colour=Cluster, label=Ortholog),
  layout="rectangular"
)
gp$data = left_join(gp$data, select(ortholog_clusters, -Cluster))
gp$data = left_join(gp$data, select(ortholog_mrca, Cluster, node), by="node")

# Rename clusters according to order in plot
new_clusters = gp$data %>%
  filter(node %in% ortholog_mrca$MRCA) %>%
  select(node, y) %>%
  rename(MRCA = node) %>%
  inner_join(ortholog_mrca) %>%
  select(Cluster, y) %>%
  distinct() %>%
  arrange(-y) %>%
  mutate(NewCluster = as.character(1:nrow(.))) %>%
  select(-y)

ortholog_clusters = ortholog_clusters %>%
  inner_join(new_clusters) %>%
  select(-Cluster) %>%
  rename(Cluster = NewCluster)

ortholog_mrca = ortholog_mrca %>%
  inner_join(new_clusters) %>%
  select(-Cluster) %>%
  rename(Cluster = NewCluster)

gp$data = gp$data %>%
  left_join(new_clusters) %>%
  select(-Cluster) %>%
  rename(Cluster = NewCluster)

gp$data = mutate(gp$data, Cluster = ifelse(is.na(Cluster), "None", Cluster))

gp$data = gp$data %>% left_join(
  ortholog_mrca %>%
    select(MRCA, Cluster) %>%
    distinct() %>%
    rename(node = MRCA, Cluster_label = Cluster)
)

gp = gp + geom_label(aes(label=Cluster_label), alpha=0.8)

gp = gp + scale_color_manual(
  values=c(
    "#d73027",
    "#f46d43",
    "#fdae61",
    "#fee090",
    "#abd9e9",
    "#74add1",
    "#4575b4",
    "#4d4d4d"
  ),
  guide = F
)
gp1 = gp

# Plot cluster properties

# Count interactions per Organism and Metabolite for each Cluster
interactions_per_cluster = interaction_comparison_concs %>%
  # Add cluster numbers
  inner_join(ortholog_clusters) %>%
  # Sum the interactions per Organism, Metabolite, and Cluster
  group_by(Organism, Metabolite, Conc, Cluster) %>%
  summarise(Interactions = sum(Interaction)) %>%
  # Add missing values as zeros for plotting
  complete(
    Organism, Metabolite, Conc, Cluster,
    fill = list(Interactions = 0)
  ) %>%
  distinct()

# Plot the interactions per Organism, Metabolite, and Cluster
gp = ggplot(
  interactions_per_cluster,
  aes(
    x=Metabolite, y=Interactions,
    group=paste(Organism, Conc),
    fill=Organism, alpha=Conc
  )
)
gp = gp + geom_col(position=position_dodge(width=0.75), width=0.75)
gp = gp + facet_grid(Cluster~., scales="free_y")
gp = gp + theme_bw()
gp = gp + theme(
  axis.text = element_text(colour="black"),
  axis.ticks = element_line(colour="black"),
  strip.background = element_blank(),
  axis.text.x = element_text(angle=60, vjust=1, hjust=1),
  strip.text = element_blank()
)
gp = gp + scale_alpha_manual(values=c(1,0.25))
gp = gp + scale_fill_manual(values=c("#9970ab","#a6dba0","#1b7837"))
gp2 = gp

# Count the number of interactions per Organism, Ortholog Category, and Cluster
interactions_per_category = interaction_comparison_concs %>%
  # Add cluster numbers
  inner_join(ortholog_clusters) %>%
  # Add eggNOG Category annotations
  inner_join(eggnog_annotations_unique) %>%
  # Sum the interactions per Organism, Category, and Cluster
  group_by(Organism, Conc, Category, Cluster) %>%
  summarise(Interactions = sum(Interaction)) %>%
  ungroup() %>%
  # Add missing values as zero
  complete(Organism, Category, Conc, Cluster, fill = list(Interactions = 0)) %>%
  distinct() %>%
  # Add eggNOG Category Description
  inner_join(select(eggnog_category_descriptions, -Description)) %>%
  # Create new Label of Category Description
  mutate(Label = paste("[", Category, "] ", Short_description, sep="")) %>%
  # Order Label by the number of Interactions
  mutate(
    Label = factor(
      Label,
      levels = (
        group_by(., Label) %>%
        summarise(Count = sum(Interactions)) %>%
        arrange(-Count, Label) %>%
        pull(Label)
      )
    )
  )

# Plot the interactions per Organism, Category, and Cluster
gp = ggplot(
  interactions_per_category,
  aes(
    x=Label, y=Interactions,
    group=paste(Organism, Conc),
    fill=Organism, alpha=Conc
  )
)
gp = gp + geom_col(position=position_dodge(width=0.75), width=0.75)
gp = gp + facet_grid(Cluster~., scales="free_y")
gp = gp + theme_bw()
gp = gp + theme(
  axis.text = element_text(colour="black"),
  axis.ticks = element_line(colour="black"),
  strip.background = element_blank(),
  axis.text.x = element_text(angle=60, vjust=1, hjust=1),
  axis.title.y = element_blank()
)
gp = gp + scale_alpha_manual(values=c(1,0.25))
gp = gp + scale_fill_manual(values=c("#9970ab","#a6dba0","#1b7837"))
gp = gp + xlab("Ortholog category")
gp3 = gp

# Make combined plot
outfile = "results/orthologs_interaction_clustering.pdf"
pdf(outfile, width=12.5, height=8.5, onefile=FALSE)
ggarrange(
  gp1,
  ggarrange(
    gp2,
    gp3,
    nrow = 1,
    ncol = 2,
    labels = c("B", ""),
    widths = c(1, 1),
    common.legend = T,
    legend="bottom",
    align="h"
  ),
  common.legend=F,
  nrow = 1,
  ncol = 2,
  labels = c("A",""),
  widths = c(1,4)
)
garbage = dev.off()

# Calculate and plot Ortholog Categories per Metabolite in total
interactions_per_category_and_metabolite = interaction_comparison_concs %>%
  # Add eggNOG Category annotations
  inner_join(eggnog_annotations_unique) %>%
  # Count number of Interactions per Organism, Metabolite, and Category
  group_by(Organism, Metabolite, Conc, Category) %>%
  summarise(Interactions = sum(Interaction)) %>%
  ungroup() %>%
  # Add missing data as zero for plotting
  complete(
    Organism, Metabolite, Conc, Category, fill = list(Interactions = 0)
  ) %>%
  distinct() %>%
  # Add eggNOG Category Description
  inner_join(select(eggnog_category_descriptions, -Description)) %>%
  # Create new Label of Category Description
  mutate(Label = paste("[", Category, "] ", Short_description, sep="")) %>%
  # Order Label by the number of Interactions
  mutate(
    Label = factor(
      Label,
      levels = (
        group_by(., Label) %>%
        summarise(Count = sum(Interactions)) %>%
        arrange(-Count, Label) %>%
        pull(Label)
      )
    )
  )

# Plot Interactions summary
gp = ggplot(
  interactions_per_category_and_metabolite,
  aes(
    x=Label, y=Interactions,
    group=paste(Organism, Conc),
    fill=Organism, alpha=Conc
  )
)
gp = gp + geom_col(position=position_dodge(width=0.75), width=0.75)
gp = gp + facet_grid(~Metabolite, scales="free_x")
gp = gp + theme_bw()
gp = gp + coord_flip()
gp = gp + theme(
  axis.text = element_text(colour="black"),
  axis.ticks = element_line(colour="black"),
  strip.background = element_blank(),
  strip.text.y = element_text(angle=0, vjust=0.5, hjust=0),
  axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)
)
gp = gp + scale_alpha_manual(values=c(1,0.25))
gp = gp + scale_fill_manual(values=c("#9970ab","#a6dba0","#1b7837"))
gp = gp + xlab("Ortholog category")

ggsave("results/metabolite_function_interactions.pdf", gp, w=16, h=7)

# Now cluster based on metabolite or ortholog for each organism
cluster_organism = function(organism, grouping, n_clusters){

  # Define cluster color palette
  available_colors = c(
    "#a50026", "#d73027", "#f46d43", "#fdae61",
    "#abd9e9", "#74add1", "#4575b4", "#313695"
  )
  custom_palette = c(
    available_colors[c(
      1:(ceiling(n_clusters/2)),
      (
        length(available_colors) - floor(n_clusters/2) + 1
      ):length(available_colors)
    )],
    "#4d4d4d"
  )

  # Prepare distance matrix
  jaccard_dist = interaction_comparison_concs_wide %>%
    filter(Organism == organism) %>%
    select(-Organism, -Metabolite, -Conc) %>%
    as.matrix()

  # Remove missing orthologs
  jaccard_dist = jaccard_dist[,colSums(jaccard_dist, na.rm=T) > 0]

  # Treat distance object differently depending on grouping
  if (grouping == "Ortholog"){
    # Transpose if clustering orthologs
    jaccard_dist = t(jaccard_dist)
  } else {
    # Add rownames if clustering metabolites
    rownames(jaccard_dist) = interaction_comparison_concs_wide %>%
      filter(Organism == organism) %>%
      mutate(Met_Conc = paste(Metabolite, Conc, sep="_")) %>%
      pull(Met_Conc)
  }

  # Calculate Jaccard distance
  jaccard_dist = vegdist(jaccard_dist, "jaccard", na.rm=T)

  # Replace missing values with maximum distance (1)
  jaccard_dist_m = as.matrix(jaccard_dist)
  jaccard_dist_m[!is.finite(jaccard_dist_m)] = 1
  jaccard_dist = as.dist(jaccard_dist_m)

  # Cluster based on distance
  clustering = hclust(jaccard_dist, method = "ward.D2")
  clusters = cutree(clustering, k=n_clusters)
  clusters = tibble(
    label = names(clusters),
    Cluster = as.character(clusters)
  )

  # Create tree from clustering
  cluster_tree = as.phylo(clustering)

  # Find most recent common ancestor for each cluster, then colour by cluster
  cluster_mrca = clusters %>%
    group_by(Cluster) %>%
    summarise(MRCA = findMRCA(cluster_tree, label)) %>%
    mutate(
      node = lapply(MRCA, function(x){getDescendants(cluster_tree, x)})
    ) %>%
    unnest(cols="node")

  # Plot clustering
  gp = ggtree(
    cluster_tree,
    aes(colour=Cluster),
    layout="rectangular"
  )
  gp$data = left_join(gp$data, select(clusters, -Cluster))
  gp$data = left_join(gp$data, select(cluster_mrca, Cluster, node), by="node")

  # Rename clusters according to order in plot
  new_clusters = gp$data %>%
    filter(node %in% cluster_mrca$MRCA) %>%
    select(node, y) %>%
    rename(MRCA = node) %>%
    inner_join(cluster_mrca) %>%
    select(Cluster, y) %>%
    distinct() %>%
    arrange(-y) %>%
    mutate(NewCluster = as.character(1:nrow(.))) %>%
    select(-y)

  clusters = clusters %>%
    inner_join(new_clusters) %>%
    select(-Cluster) %>%
    rename(Cluster = NewCluster)

  cluster_mrca = cluster_mrca %>%
    inner_join(new_clusters) %>%
    select(-Cluster) %>%
    rename(Cluster = NewCluster)

  gp$data = gp$data %>%
    left_join(new_clusters) %>%
    select(-Cluster) %>%
    rename(Cluster = NewCluster)

  gp$data = mutate(gp$data, Cluster = ifelse(is.na(Cluster), "None", Cluster))

  gp$data = gp$data %>% left_join(
    cluster_mrca %>%
      select(MRCA, Cluster) %>%
      distinct() %>%
      rename(node = MRCA, Cluster_label = Cluster)
  )

  if (grouping == "Metabolite") {
    gp = gp + geom_tiplab(aes(label=label), size=2.5)
    gp = gp + xlim(0, max(gp$data$x)*1.5)
  }

  gp = gp + geom_label(aes(label=Cluster_label), alpha=0.8)
  gp = gp + scale_color_manual(
    values=custom_palette,
    guide = F
  )
  gp_tree = gp

  # Plot PCoA
  jaccard_pcoa = pcoa(jaccard_dist)

  # Summarize interactions
  interaction_summary = interaction_comparison_concs %>%
    # Filter to Organism in question
    filter(Organism == organism) %>%
    # Rename the column that is used as grouping to "label"
    mutate(Metabolite = paste(Metabolite, Conc, sep="_")) %>%
    rename(label = grouping) %>%
    # Sum the Interactions for each label
    group_by(label) %>%
    summarise(Interactions = sum(Interaction))

  # Create plotting table
  jaccard_pcoa_plot = jaccard_pcoa$vectors[,1:2] %>%
    as_tibble() %>%
    # Add label as the column names of the distance matrix
    mutate(label = colnames(jaccard_dist_m)) %>%
    # Rename axes to something sensible
    rename(PCo1 = Axis.1, PCo2 = Axis.2) %>%
    # Add the interaction counts per label
    inner_join(interaction_summary) %>%
    # Add cluster numbers
    inner_join(clusters)

  # If grouping is Metabolite, make labels shorter
  if (grouping == "Metabolite") {
    jaccard_pcoa_plot = jaccard_pcoa_plot %>%
      mutate(
        label = str_replace(label, "_High", "-H"),
        label = str_replace(label, "_Low", "-L")
      )
  }

  # Calculate fraction of variance per PC
  library(scales)
  jaccard_pcoa_var = percent(jaccard_pcoa$values$Relative_eig, accuracy=0.1)

  # Plot it
  library(ggrepel)

  gp = ggplot(
    jaccard_pcoa_plot,
    aes(x=PCo1, y=PCo2, colour=Cluster, label=label)
  )
  gp = gp + geom_point(alpha=0.6, mapping=aes(size=Interactions))

  # Add labels only if fewer than 50
  if (nrow(jaccard_pcoa_plot) < 50){
    gp = gp + geom_text_repel(force=3, size=2,alpha=0.9)
  }

  gp = gp + scale_color_manual(
    values=custom_palette,
    guide = F
  )
  gp = gp + labs(
    x=paste("PCo1 (", jaccard_pcoa_var[1], ")", sep=""),
    y=paste("PCo2 (", jaccard_pcoa_var[2], ")", sep="")
  )
  gp = gp + theme_bw()
  gp = gp + scale_size_continuous(
    guide=guide_legend(title.position = "top"),
    breaks = trans_breaks(identity, identity, n = 4)
  )
  gp = gp + theme(
    axis.ticks = element_line(colour="black"),
    axis.text = element_text(colour="black"),
    strip.background = element_blank(),
    aspect.ratio = 1,
    legend.position = "bottom"
  )
  gp_pcoa = gp

  # Plot cluster properties

  # If the grouping is "Ortholog", one may summarise interactions per metabolite
  if (grouping == "Ortholog") {
    # Count interactions per Organism and Metabolite for each Cluster
    interactions_per_cluster = interaction_comparison_concs %>%
      filter(Organism == organism) %>%
      rename(label = Ortholog) %>%
      inner_join(clusters) %>%
      group_by(Metabolite, Conc, Cluster) %>%
      summarise(Interactions = sum(Interaction)) %>%
      # Add missing values as zero for plotting
      complete(Metabolite, Conc, Cluster, fill = list(Interactions = 0)) %>%
      distinct()

    gp = ggplot(
      interactions_per_cluster,
      aes(
        x=Metabolite, y=Interactions,
        group=Conc,
        fill=Cluster, alpha=Conc
      )
    )
    gp = gp + geom_col(position=position_dodge(width=0.75), width=0.75)
    gp = gp + facet_grid(Cluster~., scales="free_y")
    gp = gp + theme_bw()
    gp = gp + theme(
      axis.text = element_text(colour="black"),
      axis.ticks = element_line(colour="black"),
      strip.background = element_blank(),
      axis.text.x = element_text(angle=60, vjust=1, hjust=1),
      strip.text = element_blank()
    )
    gp = gp + scale_fill_manual(
      values=custom_palette,
      guide = F
    )
    gp = gp + scale_alpha_manual(values=c(1,0.25))
    gp_metabolites = gp
  }

  # If grouping is Metabolite, split cluster label
  if (grouping == "Metabolite") {
    clusters = clusters %>% separate(label, c("label", "Conc"), sep="_")
  }

  # Count interactions per ortholog category
  interactions_per_category = interaction_comparison_concs %>%
    # Select relevant Organism
    filter(Organism == organism) %>%
    # Add eggNOG Categories
    inner_join(eggnog_annotations_unique) %>%
    # Rename the column that is used as grouping to "label"
    rename(label = grouping) %>%
    # Add Cluster numbers
    inner_join(clusters) %>%
    # Sum up the Interactions for each Category and Cluster
    group_by(Category, Conc, Cluster) %>%
    summarise(Interactions = sum(Interaction)) %>%
    ungroup() %>%
    # Add missing data as zero for plotting
    complete(
      Category, Conc, Cluster, fill = list(Interactions = 0)
    ) %>%
    distinct() %>%
    # Add eggNOG Category Description
    inner_join(select(eggnog_category_descriptions, -Description)) %>%
    # Create a new Label describing the Category
    mutate(Label = paste("[", Category, "] ", Short_description, sep="")) %>%
    # Sort the Label by number of Interactions
    mutate(
      Label = factor(
        Label,
        levels = levels(interactions_per_category_and_metabolite$Label)
      )
    )

  # Plot Interactions per Category ("Label") and Cluster
  gp = ggplot(
    interactions_per_category,
    aes(x=Label, y=Interactions, group=Conc, fill=Cluster, alpha=Conc)
  )
  gp = gp + geom_col(position=position_dodge(width=0.75), width=0.75)
  gp = gp + facet_grid(Cluster~., scales="free_y")
  gp = gp + theme_bw()
  gp = gp + theme(
    axis.text = element_text(colour="black"),
    axis.ticks = element_line(colour="black"),
    strip.background = element_blank(),
    axis.text.x = element_text(angle=60, vjust=1, hjust=1),
    axis.title.y = element_blank()
  )
  gp = gp + scale_fill_manual(
    values=custom_palette,
    guide = F
  )
  gp = gp + scale_alpha_manual(values=c(1,0.25))
  gp = gp + xlab("Ortholog category")
  gp_categories = gp

  # Select plot heights based on grouping
  plot_heights = if (grouping == "Ortholog") {c(2,1)} else {c(1.5,1)}

  # Make combined plot
  gp_clustering = ggarrange(
    gp_tree,
    gp_pcoa,
    nrow = 2,
    ncol = 1,
    labels = c("A", "B"),
    heights = plot_heights
  )

  if (grouping == "Ortholog") {
    gp_composition = ggarrange(
      gp_metabolites,
      gp_categories,
      nrow = 1,
      ncol = 2,
      labels = c("C", ""),
      widths = c(1, 1),
      align="h",
      common.legend=T,
      legend="bottom"
    )
  } else {
    gp_composition = ggarrange(
      gp_categories,
      nrow = 1,
      ncol = 1,
      labels = c("C", ""),
      widths = c(1),
      common.legend=T,
      legend="bottom"
    )
  }

  # Select width of plots based on grouping
  plot_widths = if (grouping == "Ortholog") {c(1,4)} else {c(1,1.5)}
  plot_width = if (grouping == "Ortholog") {11} else {11/5*3}

  outfile = paste(
    "results/ortholog_clustering.", organism, "_by_", grouping, ".pdf", sep=""
  )
  pdf(outfile, width=plot_width, height=8.5, onefile=FALSE)
  print(
    annotate_figure(
      ggarrange(
        gp_clustering,
        gp_composition,
        common.legend=F,
        nrow = 1,
        ncol = 2,
        widths = plot_widths
      ),
      top = organism
    )
  )
  dev.off()

}

# Plot each organism by Metabolite
cluster_organism("Cupriavidus", "Metabolite", 4)
cluster_organism("Synechocystis", "Metabolite", 6)
cluster_organism("Synechococcus", "Metabolite", 6)

# Plot each organism by Ortholog
cluster_organism("Cupriavidus", "Ortholog", 7)
cluster_organism("Synechocystis", "Ortholog", 7)
cluster_organism("Synechococcus", "Ortholog", 7)
