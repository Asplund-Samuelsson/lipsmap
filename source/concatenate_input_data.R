options(width=110)
library(tidyverse)

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
lipsmap_data = bind_rows(lapply(
  1:nrow(infiles),
  function(i){
    infile = infiles[i,]
    read_csv(file.path(infile$Directory, infile$File)) %>%
      mutate(
        Organism = infile$Organism,
        Metabolite = infile$Metabolite,
        Run = infile$Run,
        Date = infile$Date
      )
  }
))

# For Organism and Metabolite with multiple Date, select latest
lipsmap_data = lipsmap_data %>%
  select(Organism, Metabolite, Date) %>%
  distinct() %>%
  group_by(Organism, Metabolite) %>%
  top_n(1, Date) %>%
  ungroup() %>%
  left_join(lipsmap_data)

# Save data to compressed archive
write_tsv(
  lipsmap_data, gzfile("data/annotated_comparison_results.2021-05-17.tab.gz")
)
