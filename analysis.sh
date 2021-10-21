#!/usr/bin/env bash

# Read command line variables
INFILE=$1
OUTDIR=$2
LOGFILE="${OUTDIR}/analysis.log"

# Create output directory if it does not exist
mkdir -p $OUTDIR

# Set up analysis log file with time and date
date > $LOGFILE

# 1. Interactions with enzymes
echo -en "\n1. Interactions with enzymes\n" >> $LOGFILE
mkdir -p ${OUTDIR}/enzymes
Rscript source/enzymes.R $INFILE ${OUTDIR}/enzymes >> $LOGFILE 2>&1

# 2. Functional categories
echo -en "\n2. Functional categories\n" >> $LOGFILE
mkdir -p ${OUTDIR}/functions
Rscript source/functions.R $INFILE ${OUTDIR}/functions >> $LOGFILE 2>&1

# 3. Comparison of orthologs
echo -en "\n3. Comparison of orthologs\n" >> $LOGFILE
mkdir -p ${OUTDIR}/orthologs
Rscript source/orthologs.R $INFILE ${OUTDIR}/orthologs >> $LOGFILE 2>&1

# 4. Comparison of KEGG modules
echo -en "\n4. Comparison of KEGG modules\n" >> $LOGFILE
mkdir -p ${OUTDIR}/modules
Rscript source/modules.R $INFILE ${OUTDIR}/modules >> $LOGFILE 2>&1
Rscript source/category_module_overlap.R >> $LOGFILE 2>&1

# 5. Phylogenetic trees
echo -en "\n5. Phylogenetic trees\n" >> $LOGFILE
mkdir -p ${OUTDIR}/phylogenetics/ko_trees
Rscript source/phylogenetics.R $INFILE ${OUTDIR}/phylogenetics >> $LOGFILE 2>&1

# 6. Supplementary tables
echo -en "\n6. Supplementary tables\n" >> $LOGFILE
mkdir -p ${OUTDIR}/tables
Rscript source/tables.R $INFILE ${OUTDIR}/tables >> $LOGFILE 2>&1

# 7. Calvin cycle and sinks
echo -en "\n7. Calvin cycle and sinks\n" >> $LOGFILE
mkdir -p ${OUTDIR}/cbb
Rscript source/cbb.R $INFILE ${OUTDIR}/cbb >> $LOGFILE 2>&1

# 8. Carbon concentration mechanisms
echo -en "\n8. Carbon concentration mechanisms\n" >> $LOGFILE
mkdir -p ${OUTDIR}/ccm
Rscript source/ccm.R $INFILE ${OUTDIR}/ccm >> $LOGFILE 2>&1

# 9. Number of detected peptides
echo -en "\n9. Number of detected peptides\n" >> $LOGFILE
mkdir -p ${OUTDIR}/peptides
Rscript source/peptides.R $INFILE ${OUTDIR}/peptides >> $LOGFILE 2>&1

# 10. Post-translational modifications
echo -en "\n10. Post-translational modifications\n" >> $LOGFILE
mkdir -p ${OUTDIR}/modifications
Rscript source/modifications.R $INFILE ${OUTDIR}/modifications >> $LOGFILE 2>&1

# 11. Quality control
echo -en "\n11. Quality control\n" >> $LOGFILE
mkdir -p ${OUTDIR}/qc
Rscript source/qc.R $INFILE ${OUTDIR}/qc >> $LOGFILE 2>&1
