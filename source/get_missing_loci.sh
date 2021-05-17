cat data/missing_locus_uniprot_IDs.txt | parallel --no-notice --jobs 16 '
  wget -qO - "https://www.uniprot.org/uniprot/{}.txt" |
  grep "^DR   KEGG" | tr ";:" "\t" | cut -f 3 | sed -e "s/^/{}\t/"
' > data/uniprot_locus.tab
