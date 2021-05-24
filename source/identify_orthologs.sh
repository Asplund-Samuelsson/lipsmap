# Download EggNOG data from UniProt
cat data/Cupriavidus.fasta data/Synechocystis.fasta data/Synechococcus.fasta |
grep ">" | cut -f 2 -d \| | parallel --no-notice --jobs 16 '
  wget -qO - "https://www.uniprot.org/uniprot/{}.txt" |
  grep "^DR   eggNOG" | tr ";" "\t" | cut -f 2  | sed -e "s/^/{}\t/"
' > data/uniprot_eggNOG.tab

# Create UniProt ID to organism mapping
ls data/*fasta | while read Fasta; do
  Organism=`echo "$Fasta" | tr "/." "\t" | cut -f 2`
  cat $Fasta | grep ">" | cut -f 2 -d \| | sed -e "s/^/${Organism}\t/"
done > data/organism_uniprot.tab
