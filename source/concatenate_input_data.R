options(width=110)
library(tidyverse)
library(foreach)
library(doMC)

# Define input path
inpath = "/hdd/emil/LiP-SMap/data/comparison_results/annotated"

# List infiles
infiles = tibble(
  Organism = list.files(inpath),
  Directory = list.files(inpath, full.names=T)
) %>%
  mutate(File = lapply(Directory, function(x){list.files(x)})) %>%
  unnest(File) %>%
  # Get Metabolite, Run, and Date labels
  mutate(
    Metabolite = str_remove(File, ".csv$") %>% str_split("_") %>% sapply("[",2),
    Run = str_remove(File, ".csv$") %>% str_split("_") %>% sapply("[",3),
    Date = str_remove(File, ".csv$") %>% str_split("_") %>% sapply("[",6)
  )

# Load all data
registerDoMC(32)
lipsmap = bind_rows(
  foreach(i=1:nrow(infiles)) %dopar% {
    infile = infiles[i,]
    read_csv(file.path(infile$Directory, infile$File)) %>%
      mutate(
        Organism = infile$Organism,
        Metabolite = infile$Metabolite,
        Run = infile$Run,
        Date = infile$Date
      )
  }
)

# Fix Metabolite, Conc and Locus column
lipsmap = lipsmap %>%
  rename(Locus = Gene_names_ordered_locus) %>%
  mutate(Conc = ifelse(str_ends(Metabolite, "-H"), 2, Conc)) %>%
  mutate(Metabolite = str_remove(Metabolite, "-[HL]$"))

# Save UniProt IDs of peptides with missing loci
write(
  filter(lipsmap, is.na(Locus)) %>% pull(UniProt_entry) %>% unique(),
  "data/missing_locus_uniprot_IDs.txt"
)

# Download missing loci
system(paste(c(
"cat data/missing_locus_uniprot_IDs.txt | parallel --no-notice --jobs 16 '",
'wget -qO - "https://www.uniprot.org/uniprot/{}.txt" |',
'grep "^DR   KEGG" | tr ";:" "\t" | cut -f 3 | sed -e "s/^/{}\t/"',
"' > data/uniprot_locus_missing.tab"
), collapse=""))

# Load missing locus IDs
uniprot_locus = read_tsv(
  "data/uniprot_locus_missing.tab",
  col_names = c("UniProt_entry", "Locus")
)

# Repair locus ID mapping (for Synechocystis)
lipsmap = bind_rows(
  lipsmap %>%
    filter(UniProt_entry %in% uniprot_locus$UniProt_entry) %>%
    select(-Locus) %>%
    left_join(uniprot_locus),
  lipsmap %>%
    filter(!(UniProt_entry %in% uniprot_locus$UniProt_entry))
)

# For Organism, Metabolite, and Conc with multiple Date, select latest
lipsmap = lipsmap %>%
  select(Organism, Metabolite, Conc, Date) %>%
  distinct() %>%
  group_by(Organism, Metabolite, Conc) %>%
  top_n(1, Date) %>%
  ungroup() %>%
  left_join(lipsmap)

# Keep only data of interest
lipsmap = lipsmap %>%
  select(
    Organism, Metabolite, Conc, Date, Peptide_ID,
    log2FC, SE, adj.pvalue, Peptide, UniProt_entry,
    Gene_ontology_GO, Pathway, Sign, Locus
  )

# Save data to compressed archive
write_tsv(
  lipsmap, gzfile("data/annotated_comparison_results.tab.gz")
)
