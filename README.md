# LiP-SMap data analysis

Analysis of [limited proteolysis small molecule mapping (LiP-SMap)](https://www.cell.com/cell/fulltext/S0092-8674(17)31448-4) data.

### Contents

1. [Data preparation](#prep)
   1. [Concatenation](#cat)
   2. [Additional annotations](#annotations)
   3. [Phylogenetic analysis](#phylogenetics)
2. [Analysis pipeline](#analysis)
   1. [Interaction with enzymes](#enzymes)
   2. [Functional categories](#functions)
   3. [Comparison of orthologs](#orthologs)
   4. [Comparison of KEGG modules](#modules)
   5. [Phylogenetic trees](#trees)
   6. [Supplementary tables](#tables)
   7. [Calvin cycle and sinks](#calvin)
   8. [Carbon concentration mechanisms](#ccm)
   9. [Number of detected peptides](#peptides)
   10. [Post-translational modifications](#mod)
3. [Author](#author)

<a name="prep"></a>
## 1. Data preparation

Data were prepared for analysis by concatenating input files, by gathering annotations for the organisms in question, and by performing a phylogenetic analysis. These data preparation steps were less prone to changes and updates compared to the analysis pipeline, and were therefore held separate.

<a name="cat"></a>
### 1.1. Concatenation

Data were concatenated from the original source:
```
Rscript source/concatenate_input_data.R ${INDIR} ${OUTFILE}
```

Missing loci were amended by the concatenation R script with the intermediate files `intermediate/missing_locus_uniprot_IDs.txt` and `intermediate/uniprot_locus_missing.tab`.

The R script finally produced concatenated LiP-SMap in `OUTFILE`.

<a name="annotations"></a>
### 1.2. Additional annotations

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

<a name="phylogenetics"></a>
### 1.3. Phylogenetic analysis

A phylogenetic analysis was performed on Calvin cycle enzymes from the corresponding KEGG module ([M00165](https://www.genome.jp/kegg-bin/show_module?M00165)), supplemented with transaldolase (K00616, K13810), triose-phosphate isomerase (K01803), and ribulose-phosphate epimerase (K01783). Sequences were downloaded from UniProt, filtered with CD-HIT, aligned with MAFFT, and used to make trees with FastTreeMP:
```
source/Calvin_cycle_phylogenetics.sh
```

Within the phylogenetic analysis script, the NCBI taxonomy was consulted to give organism labels to all proteins in the tree using the helper script `source/taxid-to-taxonomy.py`. Furthermore, UniProt sequences were acquired based on KEGG orthologs using the helper script `source/uniprot_sequences_from_KO.sh`.

<a name="analysis"></a>
## 2. Analysis pipeline

A series of R scripts were used to test hypotheses, and to generate data visualizations and tables. All scripts can be run in order by invoking the analysis shell script:

```
./analysis.sh ${INFILE} ${OUTDIR}
```

...where `INFILE` is an output file from `source/concatenate_input_data.R`, and `OUTDIR` is a desired output directory, in which each script will have its own output subdirectory.

Standard output ("`stdout`") and standard error ("`stderr`") reporting from the scripts is logged in this file:

```
${OUTDIR}/analysis.log
```

The individual steps of the analysis are described below.

<a name="enzymes"></a>
### 2.1. Interactions with enzymes

Metabolite interactions with enzymes and non-enzymes were compared using Fisher's exact test:
```
Rscript source/enzymes.R ${INFILE} ${OUTDIR}/enzymes
```

...yielding the following results:
```
Fisher_exact_test_for_enzyme_interactions.tab
```

<a name="functions"></a>
### 2.2. Functional categories

Proteins were grouped by various functional categories (EC, GO, KEGG module, pathway) and tested for enrichment of interactions:

```
Rscript source/functions.R ${INFILE} ${OUTDIR}/functions
```

...yielding the following results:
```
Fisher_exact_test_for_EC_interactions.tab
Fisher_exact_test_for_GO_interactions.tab
Fisher_exact_test_for_module_interactions.tab
Fisher_exact_test_for_pathway_interactions.tab
```

<a name="orthologs"></a>
### 2.3. Comparison of orthologs

Interaction patterns with orthologs were compared within and between organisms:
```
Rscript source/orthologs.R ${INFILE} ${OUTDIR}/orthologs
```

...producing a range of comparison plots and tables.

Overview and clustering of orthologs across the whole dataset:
```
orthologs_interaction_comparison.png
orthologs_interaction_clustering.pdf
metabolite_function_interactions.pdf
```

PCoA and PCA analysis at high and low concentration:
```
Fig.orthologs_interaction_pcoa.high.pdf
Fig.orthologs_interaction_pcoa.low.pdf
orthologs_interaction_jaccard.high.tab
orthologs_interaction_jaccard.low.tab
orthologs_interaction_pca.high.pdf
orthologs_interaction_pca.low.pdf
orthologs_interaction_pcoa.high.pdf
orthologs_interaction_pcoa.low.pdf
```

**Example:** Ortholog metabolite interaction PCoA

![alt text](data/examples/orthologs_interaction_pcoa.png "Ortholog metabolite interaction PCoA example")

Clustering of metabolites or orthologs per organism:
```
ortholog_clustering.Cupriavidus_by_Metabolite.pdf
ortholog_clustering.Cupriavidus_by_Ortholog.pdf
ortholog_clustering.Hydrogenophaga_by_Metabolite.pdf
ortholog_clustering.Hydrogenophaga_by_Ortholog.pdf
ortholog_clustering.Synechococcus_by_Metabolite.pdf
ortholog_clustering.Synechococcus_by_Ortholog.pdf
ortholog_clustering.Synechocystis_by_Metabolite.pdf
ortholog_clustering.Synechocystis_by_Ortholog.pdf
```

**Example:** Ortholog metabolite interaction clustering in _Cupriavidus_

![alt text](data/examples/ortholog_clustering.Cupriavidus_by_Metabolite.png "Ortholog metabolite interaction clustering example in Cupriavidus")

Clustered heatmap of interactions between metabolites and ortholog categories:
```
ortholog_category_heatmap.abs.pdf
ortholog_category_heatmap.norm.pdf
```

<a name="modules"></a>
### 2.4. Comparison of KEGG modules

KEGG modules are groups of enzymes constituting complete or partial pathways. Proteins were grouped by these modules and interactions were summarized and compared:
```
Rscript source/modules.R ${INFILE} ${OUTDIR}/modules
```

...yielding the following results:
```
module_interaction_summary.tab
module_interactions.pdf
```

The overlap of modules and ortholog categories was examined:
```
Rscript source/category_module_overlap.R
```

...producing the following plot:
```
results/category_module_overlap.pdf
```

**Example:** KEGG module metabolite interactions (top modules by number of interactions)

![alt text](data/examples/module_interactions.png "KEGG module metabolite interactions (top modules)")

<a name="trees"></a>
### 2.5. Phylogenetic trees

Phylogenetic trees of Calvin cycle genes were plotted using _phytools_ and _ggtree_ in R:
```
Rscript source/phylogenetics.R ${INFILE} ${OUTDIR}/phylogenetics
```

...producing the following final PDF containing visualizations of all trees, highlighting interactions with metabolites:
```
cbb_ko_trees.pdf
```

**Example:** PRK phylogenetic tree with LiP-SMap interactions

![alt text](data/examples/cbb_ko_trees.png "PRK phylogenetic tree")

<a name="tables"></a>
### 2.6. Supplementary tables

KEGG EC and module annotations were combined with eggNOG orthologs and categories to give context to metabolite-protein interactions for low and high concentration:

```
Rscript source/tables.R ${INFILE} ${OUTDIR}/tables
```

...yielding a long format supplementary table:

```
ortholog_ec_module_interactions.tab.gz
ortholog_ec_module_interactions.xlsx
```

<a name="calvin"></a>
### 2.7. Calvin cycle and sinks

Calvin (CBB) cycle enzymes and related sink (or "drain") enzymes were investigated for interactions with the tested metabolites:

```
Rscript source/cbb.R ${INFILE} ${OUTDIR}/cbb
```

...by creating a plot of all interactions:

```
cbb.pdf
```

...with enzymes defined by EC number annotations from KEGG in the files `data/cbb_enzymes.tab` and `data/cbb_drains.tab`.

<a name="ccm"></a>
### 2.8. Carbon concentration mechanisms

Carbon concentration mechanism (CCM) interactions in _Synechocystis_ were plotted using this script:

```
Rscript source/ccm.R ${INFILE} ${OUTDIR}/ccm
```

...generating this output:

```
ccm.pdf
```

...based on CCM and regulatory genes listed in `data/CCM_regulatory_proteins.csv`.

<a name="peptides"></a>
### 2.9. Number of detected peptides

The number of detected peptides for proteins where none (Interaction FALSE) or at least one peptide (Interaction TRUE) showed significant interaction with a metabolite was plotted for each organism and concentration:

```
Rscript source/peptides.R ${INFILE} ${OUTDIR}/peptides
```

...yielding this graphic:

```
peptides_per_protein.pdf
```

**Example:** Detected peptides per protein compared to occurrence of interaction in _Synechocystis_

![alt text](data/examples/peptides.png "Detected peptides and interaction in Synechocystis")

Furthermore, the script counted the number of detected peptides in each experiment:

```
peptides_per_experiment.pdf
```

<a name="mod"></a>
### 2.10. Post-translational modifications

Interactions between _Synechocystis_ proteins and acetyl-CoA or GAP were compared to published [acetylome](https://pubs.acs.org/doi/10.1021/pr501275a) and [propionylome](https://www.mdpi.com/1422-0067/20/19/4792/html) data using Fisher's exact test followed by plotting:

```
Rscript source/modifications.R ${INFILE} ${OUTDIR}/modifications
```

...with this plot as the result:

```
modifications.pdf
```

**Example:** Overlap between post-translational modifications and interactions with acetyl-CoA or GAP in _Synechocystis_ (dashed lines indicate insignificant enrichment)

![alt text](data/examples/modifications.png "Overlap between interaction and post-translational modifications")

<a name="author"></a>
## 3. Author
Johannes Asplund-Samuelsson, KTH (johannes.asplund.samuelsson@scilifelab.se)
