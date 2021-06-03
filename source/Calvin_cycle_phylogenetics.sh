# Get all KEGG Orthologs for the Calvin cycle
wget -qO - http://rest.kegg.jp/link/ko/M00165 | cut -f 3 -d : | \
sort | uniq > data/cbb_ko.txt

# Get descriptions for Calvin cycle orthologs
wget -qO - http://rest.kegg.jp/list/`cat data/cbb_ko.txt | tr "\n" "+"` \
> data/cbb_ko_description.tab

# Download UniProt sequences for each Calvin cycle KO
cat data/cbb_ko.txt | while read KO; do
  source/uniprot_sequences_from_KO.sh -k ${KO} -o data/cbb_ko/${KO}.fasta
done

# Downsample sequences based on percent identity
cd-hit -T 0 -c 1 -n 3 -i data/cbb_ko/K00615.fasta \
-o intermediate/K00615.cdhit.fasta

ls data/cbb_ko/ | while read Fasta; do
  # Cluster at decreasing identity thresholds
  for c in 1 0.95 0.9 0.85 0.8 0.75 0.7 0.65 0.6 0.55 0.5; do
    # Select word length based on threshold (why is this not automatic..?)
    if (( $(echo "$c >= 0.7" | bc -l) )); then
      n="5"
    elif (( $(echo "$c >= 0.6" | bc -l) )); then
      n="4"
    else
      n="3"
    fi
    # Perform cd-hit clustering
    cd-hit -T 0 -c $c -n $n -i data/cbb_ko/${Fasta} \
    -o intermediate/cdhit.${Fasta}
    # If there are fewer than 1000 sequences, stop clustering
    if [[ `grep -c ">" intermediate/cdhit.${Fasta}` -lt 1000 ]]; then
      mv intermediate/cdhit.${Fasta} \
        data/cbb_ko_clst/cdhit_c${c}_n${n}.${Fasta}
      break
    fi
  done
done

# Perform alignment and tree calculation
cat data/cbb_ko.txt | while read KO; do

  # Get missing sequence IDs after clustering
  (
    # Determine existing organism sequences for KO
    (
      cat data/organism_uniprot.tab | cut -f 2
      grep ">" data/cbb_ko/${KO}.fasta | cut -f 2 -d \|
    ) | sort | uniq -d
    # Determine organism sequences for KO after clustering
    (
      cat data/organism_uniprot.tab | cut -f 2
      grep ">" data/cbb_ko_clst/*.${KO}.fasta | cut -f 2 -d \|
    ) | sort | uniq -d
    # Determine missing sequences after clustering
  ) | sort | uniq -u | while read ID; do
    grep -P "\|${ID}\|" data/cbb_ko/${KO}.fasta | cut -f 1 -d \  | tr -d ">"
  done > intermediate/missing_${KO}.txt

  # Align KO sequences
  (
    seqmagick convert --include-from-file intermediate/missing_${KO}.txt \
      --output-format fasta data/cbb_ko/${KO}.fasta -
    cat data/cbb_ko_clst/*.${KO}.fasta
  ) | mafft --thread 24 - > intermediate/ko_ali/${KO}.ali.fasta

  # Create tree for KO sequences
  fasttreeMP intermediate/ko_ali/${KO}.ali.fasta > results/ko_trees/${KO}.tree

done

# Download NCBI taxonomy
cd data/ncbi/taxonomy
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -zxvf taxdump.tar.gz names.dmp nodes.dmp
rm taxdump.tar.gz
cd ../../..

# Get organism information for members of trees
cat data/cbb_ko.txt | while read KO; do

  # Obtain taxonomy IDs and taxonomy information
  (grep ">" intermediate/ko_ali/${KO}.ali.fasta | tr " " "\n" | \
  grep -P "^>|^OX=" | sed -e 's/^OX=//' | tr ">" "&" | tr "\n" "\t" | \
  tr "&" "\n" | sed -e 's/\t$//' | grep -v "^$") > \
  intermediate/ko_tax/${KO}_taxids.tab

  # Getting full taxonomy information for the taxonomy IDs
  source/taxid-to-taxonomy.py \
  -i intermediate/ko_tax/${KO}_taxids.tab \
  -n data/ncbi/taxonomy/names.dmp \
  -d data/ncbi/taxonomy/nodes.dmp \
  -o intermediate/ko_tax/${KO}_taxonomy.tab

done
