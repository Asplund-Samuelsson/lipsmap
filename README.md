# LiP-SMap data analysis

Analysis of [limited proteolysis small molecule mapping (LiP-SMap)](https://www.cell.com/cell/fulltext/S0092-8674(17)31448-4) data.

### Contents

**1. [Data preparation](#prep)**

**2. [Interaction with enzymes](#enzymes)**

**3. [Functional categories](#functions)**

**4. [Comparison of orthologs](#orthologs)**

**5. [Comparison of KEGG modules](#modules)**

**6. [Phylogenetic trees](#phylogenetics)**

**7. [Supplementary tables](#tables)**

**8. [Calvin cycle and sinks](#calvin)**

**9. [Carbon concentration mechanisms](#ccm)**

**10. [Number of detected peptides](#peptides)**

**A. [Author](#author)**

<a name="prep"></a>
## 1. Data preparation

### Concatenation

Data were concatenated from the original source:
```
source/concatenate_input_data.R
```

Missing loci were amended:
```
data/missing_locus_uniprot_IDs.txt
data/uniprot_locus_missing.tab
```

The R script finally produced concatenated LiP-SMap results:
```
data/annotated_comparison_results.tab.gz
```

### Additional annotations

#### KEGG annotations

Additional annotations were downloaded from KEGG:
```
source/get_ko_and_modules.sh
```

...providing these annotation files:
```
data/KEGGgene_uniprot_organism.tab
data/KEGGgene_KO_organism.tab
data/KEGGgene_module_organism.tab
data/KEGGgene_EC_organism.tab
data/EC_description.tab
data/module_description.tab
data/module_compound.tab
data/module_reaction.tab
data/reaction_compound.tab
```

#### eggNOG ortholog annotations

Orthologs in the organisms were identified via UniProt and the eggNOG labels:
```
source/identify_orthologs.sh
```

...resulting in these data files:
```
data/uniprot_eggNOG.tab
data/organism_uniprot.tab
```

The second file lists all UniProt sequences in the organisms, so that those without an ortholog could be accounted for as well.

Annotations for the eggNOG labels were acquired:
```
source/get_eggNOG_categories.sh
```

...and filtered with R:
```
source/filter_eggNOG_annotations.R
```

...thus producing the final eggNOG ortholog annotations:
```
data/eggNOG_annotations.tab
```

#### Gene identifiers

Gene identifiers were downloaded for the LiP-SMap organism proteins:
```
source/get_uniprot_gene_identifiers.sh
```

...producing these lists:
```
data/uniprot_gene.tab
data/uniprot_locus_complete.tab
```

### Phylogenetic analysis

A phylogenetic analysis was performed on Calvin cycle enzymes from the corresponding KEGG module ([M00165](https://www.genome.jp/kegg-bin/show_module?M00165)), supplemented with transaldolase (K00616, K13810), triose-phosphate isomerase (K01803), and ribulose-phosphate epimerase (K01783). Sequences were downloaded from UniProt, filtered with CD-HIT, aligned with MAFFT, and used to make trees with FastTreeMP:
```
source/Calvin_cycle_phylogenetics.sh
```

The NCBI taxonomy was consulted to give organism labels to all proteins in the tree using the helper script `source/taxid-to-taxonomy.py`.

UniProt sequences were acquired based on KEGG orthologs using the helper script `source/uniprot_sequences_from_KO.sh`.

<a name="enzymes"></a>
## 2. Interactions with enzymes

Metabolite interactions with enzymes and non-enzymes were compared using Fisher's exact test:
```
source/enzymes.R
```

...yielding the following results:
```
results/Fisher_exact_test_for_enzyme_interactions.tab
```

<a name="functions"></a>
## 3. Functional categories

Proteins were grouped by various functional categories (EC, GO, KEGG module, pathway) and tested for enrichment of interactions:

```
source/functions.R
```

...yielding the following results:
```
results/Fisher_exact_test_for_EC_interactions.tab
results/Fisher_exact_test_for_GO_interactions.tab
results/Fisher_exact_test_for_module_interactions.tab
results/Fisher_exact_test_for_pathway_interactions.tab
```

<a name="orthologs"></a>
## 4. Comparison of orthologs

Interaction patterns with orthologs were compared within and between organisms:
```
source/orthologs.R
```

...producing a range of comparison plots and tables.

Overview and clustering of orthologs across the whole dataset:
```
results/orthologs_interaction_comparison.png
results/orthologs_interaction_clustering.pdf
results/metabolite_function_interactions.pdf
```

PCoA and PCA analysis at high and low concentration:
```
results/Fig.orthologs_interaction_pcoa.high.pdf
results/Fig.orthologs_interaction_pcoa.low.pdf
results/orthologs_interaction_jaccard.high.tab
results/orthologs_interaction_jaccard.low.tab
results/orthologs_interaction_pca.high.pdf
results/orthologs_interaction_pca.low.pdf
results/orthologs_interaction_pcoa.high.pdf
results/orthologs_interaction_pcoa.low.pdf
```

**Example:** Ortholog metabolite interaction PCoA

![alt text](data/examples/orthologs_interaction_pcoa.png "Ortholog metabolite interaction PCoA example")

Clustering of metabolites or orthologs per organism:
```
results/ortholog_clustering.Cupriavidus_by_Metabolite.pdf
results/ortholog_clustering.Cupriavidus_by_Ortholog.pdf
results/ortholog_clustering.Hydrogenophaga_by_Metabolite.pdf
results/ortholog_clustering.Hydrogenophaga_by_Ortholog.pdf
results/ortholog_clustering.Synechococcus_by_Metabolite.pdf
results/ortholog_clustering.Synechococcus_by_Ortholog.pdf
results/ortholog_clustering.Synechocystis_by_Metabolite.pdf
results/ortholog_clustering.Synechocystis_by_Ortholog.pdf
```

**Example:** Ortholog metabolite interaction clustering in _Cupriavidus_

![alt text](data/examples/ortholog_clustering.Cupriavidus_by_Metabolite.png "Ortholog metabolite interaction clustering example in Cupriavidus")

Clustered heatmap of interactions between metabolites and ortholog categories:
```
results/ortholog_category_heatmap.abs.pdf
results/ortholog_category_heatmap.norm.pdf
```

<a name="modules"></a>
## 5. Comparison of KEGG modules

KEGG modules are groups of enzymes constituting complete or partial pathways. Proteins were grouped by these modules and interactions were summarized and compared:
```
source/modules.R
```

...yielding the following results:
```
results/module_interaction_summary.tab
results/module_interactions.pdf
```

The overlap of modules and ortholog categories was examined:
```
source/category_module_overlap.R
```

...producing the following plot:
```
results/category_module_overlap.pdf
```

**Example:** KEGG module metabolite interactions (top modules by number of interactions)

![alt text](data/examples/module_interactions.png "KEGG module metabolite interactions (top modules)")

<a name="phylogenetics"></a>
## 6. Phylogenetic trees

Phylogenetic trees of Calvin cycle genes were plotted using _phytools_ and _ggtree_ in R:
```
source/phylogenetics.R
```

...producing the following final PDF containing visualizations of all trees, highlighting interactions with metabolites:
```
results/cbb_ko_trees.pdf
```

**Example:** PRK phylogenetic tree with LiP-SMap interactions

![alt text](data/examples/cbb_ko_trees.png "PRK phylogenetic tree")

<a name="tables"></a>
## 7. Supplementary tables

KEGG EC and module annotations were combined with eggNOG orthologs and categories to give context to metabolite-protein interactions for low and high concentration:

```
source/tables.R
```

...yielding a long format supplementary table:

```
results/ortholog_ec_module_interactions.tab.gz
```

<a name="calvin"></a>
## 8. Calvin cycle and sinks

Calvin (CBB) cycle enzymes and related sink (or "drain") enzymes were investigated for interactions with the tested metabolites:

```
source/cbb_drains.R
```

...by creating a plot of all interactions:

```
results/cbb_drains.pdf
```

...with these enzymes, defined by EC number annotations from KEGG:

```
data/cbb_enzymes.tab
data/cbb_drains.tab
```

<a name="ccm"></a>
## 9. Carbon concentration mechanisms

Carbon concentration mechanism (CCM) interactions in _Synechocystis_ were plotted using this script:

```
source/ccm.R
```

...generating this output:

```
results/ccm.pdf
```

...based on this list of CCM and regulatory genes:

```
data/CCM_regulatory_proteins.csv
```

<a name="peptides"></a>
## 10. Number of detected peptides

The number of detected peptides for proteins where none (Interaction FALSE) or at least one peptide (Interaction TRUE) showed significant interaction with a metabolite was plotted for each organism and concentration:

```
source/peptides.R
```

...yielding this graphic:

```
results/peptides.pdf
```

**Example:** Detected peptides per protein compared to occurrence of interaction in _Synechocystis_

![alt text](data/examples/peptides.png "Detected peptides and interaction in Synechocystis")

<a name="author"></a>
## A. Author
Johannes Asplund-Samuelsson, KTH (johannes.asplund.samuelsson@scilifelab.se)
