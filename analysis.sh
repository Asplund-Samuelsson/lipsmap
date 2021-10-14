#!/usr/bin/env bash

# Set up analysis log file with time and date
date > analysis.log

# 1. Interactions with enzymes
echo -en "\n1. Interactions with enzymes\n" >> analysis.log
Rscript source/enzymes.R >> analysis.log 2>&1

# 2. Functional categories
echo -en "\n2. Functional categories\n" >> analysis.log
Rscript source/functions.R >> analysis.log 2>&1

# 3. Comparison of orthologs
echo -en "\n3. Comparison of orthologs\n" >> analysis.log
Rscript source/orthologs.R >> analysis.log 2>&1

# 4. Comparison of KEGG modules
echo -en "\n4. Comparison of KEGG modules\n" >> analysis.log
Rscript source/modules.R >> analysis.log 2>&1
Rscript source/category_module_overlap.R >> analysis.log 2>&1

# 5. Phylogenetic trees
echo -en "\n5. Phylogenetic trees\n" >> analysis.log
Rscript source/phylogenetics.R >> analysis.log 2>&1

# 6. Supplementary tables
echo -en "\n6. Supplementary tables\n" >> analysis.log
Rscript source/tables.R >> analysis.log 2>&1

# 7. Calvin cycle and sinks
echo -en "\n7. Calvin cycle and sinks\n" >> analysis.log
Rscript source/cbb_drains.R >> analysis.log 2>&1

# 8. Carbon concentration mechanisms
echo -en "\n8. Carbon concentration mechanisms\n" >> analysis.log
Rscript source/ccm.R >> analysis.log 2>&1

# 9. Number of detected peptides
echo -en "\n9. Number of detected peptides\n" >> analysis.log
Rscript source/peptides.R >> analysis.log 2>&1

# 10. Post-translational modifications
echo -en "\n10. Post-translational modifications\n" >> analysis.log
Rscript source/modifications.R >> analysis.log 2>&1
