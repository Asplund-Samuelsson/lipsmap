# Download EggNOG data from UniProt
cat data/*.fasta |
grep ">" | cut -f 2 -d \| | parallel --no-notice --jobs 16 '
  wget -qO - "https://www.uniprot.org/uniprot/{}.txt" |
  grep "^DR   eggNOG" | tr ";" "\t" | cut -f 2  | sed -e "s/^/{}\t/" | tr -d " "
' > data/uniprot_eggNOG.tab

# Add Hydrogenophaga orthologs manually (not available in UniProt)
cat data/Hydrogenophaga.MM_2rpcpd43.emapper.annotations.tab | grep -v "^#" |
cut -f 1,5 | tr "|@" "\t" | cut -f 2,4 >> data/uniprot_eggNOG.tab

# Create UniProt ID to organism mapping
ls data/*fasta | while read Fasta; do
  Organism=`echo "$Fasta" | tr "/." "\t" | cut -f 2`
  cat $Fasta | grep ">" | cut -f 2 -d \| | sed -e "s/^/${Organism}\t/"
done > data/organism_uniprot.tab
