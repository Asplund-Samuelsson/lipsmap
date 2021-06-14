options(width=110)
library(tidyverse)

# Define organisms and colors
organisms = c("Hydrogenophaga", "Cupriavidus", "Synechococcus", "Synechocystis")
organcols = c("#762a83", "#9970ab","#5aae61","#1b7837")

# Load data
uniprot_genename = read_tsv(
  "data/uniprot_gene.tab", col_names=c("UniProt_entry", "Gene_name")
)

uniprot_locus = read_tsv(
  "data/uniprot_locus.tab", col_names=c("UniProt_entry", "Locus")
)

org_uniprot = read_tsv(
  "data/organism_uniprot.tab",
  col_names = c("Organism", "UniProt_entry")
)

orthologs = scan("data/cbb_ko.txt", character())

ko_description = read_tsv(
  "data/cbb_ko_description.tab", col_names=c("Ortholog", "Description")
) %>% mutate(Ortholog = str_remove(Ortholog, "ko:"))

ko_name = read_tsv("data/cbb_ko_name.tab", col_names=c("Ortholog", "Name"))

taxonomy = bind_rows(
  lapply(
    orthologs,
    function(ko){
      read_tsv(paste("intermediate/ko_tax/", ko, "_taxonomy.tab", sep="")) %>%
        mutate(Ortholog = ko)
    }
  )
)

lipsmap = read_tsv(
  "data/annotated_comparison_results.tab.gz",
  col_types = cols(Protein.names = col_character())
)

interactions = lipsmap %>%
  # Determine metabolite-interacting proteins; at least one significant peptide
  group_by(Organism, Metabolite, Conc, UniProt_entry) %>%
  summarise(Interaction = "sign" %in% Sign) %>%
  # Clarify concentration description
  mutate(Conc = ifelse(Conc == 1, "Low", "High"))

# Remove orthologs that have no representatives in the organisms of interest
taxonomy = taxonomy %>%
  mutate(UniProt_entry = sapply(str_split(identifier, "\\|"), "[[", 2)) %>%
  group_by(Ortholog) %>%
  filter(length(intersect(UniProt_entry, org_uniprot$UniProt_entry)) > 0)

# Determine interactions to show in tree
liplab = interactions %>%
  semi_join(taxonomy) %>%
  filter(Interaction) %>%
  ungroup() %>%
  select(UniProt_entry, Metabolite) %>%
  distinct() %>%
  group_by(UniProt_entry) %>%
  summarise(Interaction = str_wrap(paste(Metabolite, collapse=", "), 25))

# Load trees
library(phytools)
library(ggtree)

tree_orthologs = taxonomy %>% pull(Ortholog) %>% unique()

trees = lapply(
  tree_orthologs,
  function(ko){read.tree(paste("results/ko_trees/", ko, ".tree", sep=""))}
)

# Determine unified color scheme
fcut = 0.02

tax_colors = taxonomy %>%
  select(Ortholog, identifier, group, superkingdom, kingdom, phylum) %>%
  mutate(
    Group = case_when(
      superkingdom == "Eukaryota" & !is.na(kingdom) ~ kingdom,
      superkingdom == "Viruses" ~ superkingdom,
      is.na(group) & !is.na(kingdom) ~ kingdom,
      is.na(group) ~ superkingdom,
      T ~ group
    )
  ) %>%
  group_by(Group) %>%
  mutate(Fraction = length(identifier) / nrow(.)) %>%
  ungroup() %>%
  mutate(
    Group = ifelse(
      Fraction < fcut,
      case_when(
        superkingdom == "Eukaryota" & is.na(kingdom) ~ "Other Eukaryota",
        superkingdom == "Eukaryota" ~ paste("Other", kingdom),
        phylum == "Proteobacteria" ~ "Other Proteobacteria",
        superkingdom == "Bacteria" ~ "Other Bacteria",
        superkingdom == "Archaea" ~ "Other Archaea",
        T ~ "Other"
      ),
      Group
    )
  ) %>%
  group_by(Group) %>%
  mutate(Fraction = length(identifier) / nrow(.)) %>%
  ungroup()

# Set up colors
col_table = tibble(
  Group = c(
    "Euryarchaeota", "Other Archaea",
    "Alphaproteobacteria", "Gammaproteobacteria",
    "Firmicutes", "Other Bacteria",
    "Actinobacteria", "Betaproteobacteria",
    "Cyanobacteria", "Bacteroidetes",
    "Deltaproteobacteria", "Other Proteobacteria",
    "Viridiplantae", "Fungi",
    "Other Eukaryota", "Metazoa"
  ),
  Color= c(
    "#762a83", "#9970ab",
    "#543005", "#8c510a",
    "#92c5de", "#4393c3",
    "#2166ac", "#bf812d",
    "#01665e", "#053061",
    "#dfc27d", "#f6e8c3",
    "#4d9221", "#de77ae",
    "#c51b7d", "#8e0152"
  )
)

#    Group                Color
#    <chr>                <chr>
#  1 Euryarchaeota        #762a83
#  2 Other Archaea        #9970ab
#  3 Alphaproteobacteria  #543005
#  4 Gammaproteobacteria  #8c510a
#  5 Firmicutes           #92c5de
#  6 Other Bacteria       #4393c3
#  7 Actinobacteria       #2166ac
#  8 Betaproteobacteria   #bf812d
#  9 Cyanobacteria        #01665e
# 10 Bacteroidetes        #053061
# 11 Deltaproteobacteria  #dfc27d
# 12 Other Proteobacteria #f6e8c3
# 13 Viridiplantae        #4d9221
# 14 Fungi                #de77ae
# 15 Other Eukaryota      #c51b7d
# 16 Metazoa              #8e0152

tax_colors = tax_colors %>%
  left_join(col_table) %>%
  select(identifier, Group, Color) %>%
  distinct()

# Prepare organism shape scale
org_shape = tibble(Organism = organisms, Shape = c(8, 15, 16, 17))

# Prepare organism color scale
org_color = tibble(Organism = organisms, Color = organcols)

library(ggpubr)
library(ggrepel)

# Set up tree plotting function
plot_tree = function (ko) {

  # Pick tree
  ko_tree = trees[[which(tree_orthologs == ko)]]

  # Midpoint root the tree
  ko_tree = midpoint.root(ko_tree)

  # Start building circular tree plot
  gp = ggtree(ko_tree, layout="fan", size=0.2)

  # Filter taxonomic colours to organisms in Rubisco tree
  tax_colors_tree = gp$data %>%
    filter(isTip) %>%
    select(identifier = label) %>%
    left_join(select(tax_colors, identifier, Color))

  # Create dataframe with group association for heatmap
  tax_colors_tree = data.frame(
    row.names = tax_colors_tree$identifier,
    Organism = tax_colors_tree$Color,
    stringsAsFactors=F
  )

  gp = gheatmap(
    gp, tax_colors_tree,
    offset = 1, width = 0.05, colnames = F, color=NA
  )
  gp = gp + scale_fill_identity()

  # Add scale bar
  gp = gp + geom_treescale(
    x=max(gp$data$x)*1.04, y=0,
    fontsize=3, offset=3, linesize=0.2
  )

  # Add gene name and locus
  gp$data = gp$data %>%
    separate(label, c(NA, "UniProt_entry", NA), "\\|", remove=F) %>%
    left_join(uniprot_genename) %>%
    left_join(uniprot_locus) %>%
    left_join(org_uniprot) %>%
    left_join(liplab) %>%
    mutate(
      Gene = ifelse(
        is.na(Locus), NA,
        paste(
          UniProt_entry, " ", ifelse(is.na(Gene_name), "Unnamed", Gene_name),
          " (", Locus, ")\n", Interaction,
          sep=""
        )
      ),
      Organism = factor(Organism, levels=organisms)
    ) %>%
    left_join(org_shape) %>%
    left_join(org_color)

  # Add tip label and organism shape
  gp = gp + geom_tippoint(
    mapping=aes(shape=Shape, color=Color)
  )
  gp = gp + scale_shape_identity()
  gp = gp + geom_label_repel(
    mapping=aes(label=Gene, color=Color), size=2.5, fill="#ffffffcc",
    label.size=0, max.overlaps=Inf
  )
  gp = gp + scale_color_identity()

  # Remove white space around plot
  gp = gp + theme(
    plot.margin=unit(c(-3,-3,-3,-3),"cm"),
    plot.background = element_rect(fill = "transparent",colour = NA),
    # Remove color legend
    legend.position = "none"
  )

  # Save plot in new object
  gp_tree = gp

  # Create legend
  gp = ggplot(col_table, aes(y=Group, x=0, fill=Color))
  gp = gp + geom_tile()
  gp = gp + theme_bw()
  gp = gp + scale_fill_identity()
  gp = gp + theme(
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    plot.title.position = "plot"
  )
  gp = gp + coord_equal()
  gp = gp + scale_y_discrete(limits=rev)

  # Save legend in new object
  gp_legend = gp

  # Create shape legend
  gp = ggplot(
    inner_join(org_shape, org_color),
    aes(y=factor(Organism, levels=organisms), x=0, shape=Shape, color=Color)
  )
  gp = gp + geom_point()
  gp = gp + theme_bw()
  gp = gp + scale_shape_identity()
  gp = gp + scale_color_identity()
  gp = gp + theme(
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    plot.title.position = "plot",
    legend.position = "none"
  )
  gp = gp + scale_y_discrete(limits=rev)
  gp = gp + ggtitle("Organism")

  # Save shape legend in new object
  gp_shape = gp

  # Save blank plot
  gp_blank = ggplot() + theme_bw() + theme(panel.border=element_blank())

  # Write to outfile
  pdf(
    paste(
      "results/ko_trees/",
      filter(ko_name, Ortholog == ko)$Name,
      "_", ko,
      ".pdf", sep=""
    ),
    w=12, h=10, onefile=F
  )
  print(
    annotate_figure(
      ggarrange(
        gp_tree,
        ggarrange(
          gp_blank, gp_shape, gp_legend, gp_blank,
          nrow = 4, ncol = 1,
          heights=c(1,0.7,2.5,1)
        ),
        nrow = 1, ncol = 2,
        widths = c(4, 1)
      ),
      top=paste(
        filter(ko_name, Ortholog == ko)$Name,
        paste(ko, filter(ko_description, Ortholog == ko)$Description),
        sep="\n"
      )
    )
  )
  junk = dev.off()

}

library(foreach)
library(doMC)
registerDoMC(length(tree_orthologs))

junk = foreach(ko=tree_orthologs) %dopar% {plot_tree(ko)}

# Create one PDF
system(
  paste(
    "pdfunite",
    paste(
      sort(sapply(
        tree_orthologs,
        function(x){
          grep(
            x,
            grep("pdf", list.files("results/ko_trees", full.names=T), value=T),
            value=T
          )
        }
      )),
      collapse=" "
    ),
    "results/cbb_ko_trees.pdf",
    sep=" "
  )
)
