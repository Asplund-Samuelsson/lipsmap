# Download eggNOG annotations
wget -qO intermediate/eggNOG.tab \
http://eggnog5.embl.de/download/eggnog_5.0/e5.og_annotations.tsv

wget -qO intermediate/eggNOG_4.1.tab.gz \
http://eggnog5.embl.de/download/eggnog_4.1/all_OG_annotations.tsv.gz

# Filter annotations
Rscript source/filter_eggNOG_annotations.R

# Remove intermediates
rm intermediate/eggNOG.tab intermediate/eggNOG_4.1.tab.gz
