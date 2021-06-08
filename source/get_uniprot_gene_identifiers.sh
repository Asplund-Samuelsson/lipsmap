# Download gene data from UniProt
cat data/Cupriavidus.fasta data/Synechocystis.fasta data/Synechococcus.fasta |
grep ">" | cut -f 2 -d \| | parallel --no-notice --jobs 16 '
  wget -qO - "https://www.uniprot.org/uniprot/{}.txt" |
  grep "^GN   Name" | tr " " "=" | cut -f 5 -d = |
  tr -d ";" | sed -e "s/^/{}\t/"
' > data/uniprot_gene.tab

cat data/Cupriavidus.fasta data/Synechocystis.fasta data/Synechococcus.fasta |
grep ">" | cut -f 2 -d \| | parallel --no-notice --jobs 16 '
  wget -qO - "https://www.uniprot.org/uniprot/{}.txt" |
  grep "^DR   KEGG" | tr ":;" "\t" | cut -f 3 | sed -e "s/^/{}\t/"
' > data/uniprot_locus.tab
