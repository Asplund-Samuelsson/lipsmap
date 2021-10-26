options(width=110)
library(tidyverse)
library(foreach)
library(doMC)

# Define input directory
inpath = "/hdd/emil/LiP-SMap/data/comparison_results/annotated/"

# Define list of orthologs
orthologs = read_tsv(
  "data/uniprot_eggNOG.tab",
  col_names = c("UniProt_entry", "Ortholog")
) %>% inner_join(
  read_tsv(
    "data/organism_uniprot.tab",
    col_names = c("Organism", "UniProt_entry")
  )
)

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
        Date = infile$Date,
        Infile = infile$File
      )
  }
)

# Pick a random subset of 250 proteins in each organism, of same ortholog group
proteins = orthologs %>%
  # Make sure proteins were detected
  semi_join(lipsmap) %>%
  # For each ortholog, make sure at least three organisms had it
  group_by(Ortholog) %>%
  filter(sum(unique(orthologs$Organism) %in% Organism) > 2) %>%
  # Sample 250 orthologs fitting the criteria above
  ungroup() %>%
  filter(Ortholog %in% sample(.$Ortholog, 250)) %>%
  # Sample 250 proteins per organism
  group_by(Organism) %>%
  sample_n(250) %>%
  select(-Ortholog)

# Reduce to wanted metabolites and proteins
red_lipsmap = lipsmap %>%
  filter(
    str_remove(Metabolite, "-[HL]") %in%
    c("3PGA","AcCoA","ATP","Cit","FBP","G6P","GAP","Glyx","RuBP","Ru5P")
  ) %>%
  inner_join(proteins)

# Save reduced files in the output directory
outfiles = distinct(select(red_lipsmap, Organism, Infile))
for (x in 1:nrow(outfiles)) {
  outfile = outfiles[x,]
  write_csv(
    red_lipsmap %>%
      inner_join(outfile) %>%
      select(-Organism, -Metabolite, -Run, -Date, -Infile),
    file.path("data", "tutorial", "lipsmap", outfile$Organism, outfile$Infile)
  )
}
