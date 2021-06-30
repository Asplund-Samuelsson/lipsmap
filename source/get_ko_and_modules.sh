# Get UniProt to KEGG gene ID conversion
(
echo -e "KEGGgene\tUniProt\tOrganism"
wget -qO - http://rest.kegg.jp/conv/uniprot/reh | sed -e 's/$/\tCupriavidus/'
wget -qO - http://rest.kegg.jp/conv/uniprot/syf | sed -e 's/$/\tSynechococcus/'
wget -qO - http://rest.kegg.jp/conv/uniprot/syn | sed -e 's/$/\tSynechocystis/'
wget -qO - http://rest.kegg.jp/conv/uniprot/hpse | sed -e 's/$/\tHydrogenophaga/'
) > data/KEGGgene_uniprot_organism.tab

# Get KEGG gene ID to KO conversion
(
echo -e "KEGGgene\tKO\tOrganism"
wget -qO - http://rest.kegg.jp/link/ko/reh | sed -e 's/$/\tCupriavidus/'
wget -qO - http://rest.kegg.jp/link/ko/syf | sed -e 's/$/\tSynechococcus/'
wget -qO - http://rest.kegg.jp/link/ko/syn | sed -e 's/$/\tSynechocystis/'
wget -qO - http://rest.kegg.jp/link/ko/hpse | sed -e 's/$/\tHydrogenophaga/'
) > data/KEGGgene_KO_organism.tab

# KEGG gene ID to Module conversion
(
echo -e "KEGGgene\tModule\tOrganism"
wget -qO - http://rest.kegg.jp/link/module/reh | sed -e 's/$/\tCupriavidus/'
wget -qO - http://rest.kegg.jp/link/module/syf | sed -e 's/$/\tSynechococcus/'
wget -qO - http://rest.kegg.jp/link/module/syn | sed -e 's/$/\tSynechocystis/'
wget -qO - http://rest.kegg.jp/link/module/hpse | sed -e 's/$/\tHydrogenophaga/'
) > data/KEGGgene_module_organism.tab

# KEGG gene ID to EC conversion
(
echo -e "KEGGgene\tEC\tOrganism"
wget -qO - http://rest.kegg.jp/link/ec/reh | sed -e 's/$/\tCupriavidus/'
wget -qO - http://rest.kegg.jp/link/ec/syf | sed -e 's/$/\tSynechococcus/'
wget -qO - http://rest.kegg.jp/link/ec/syn | sed -e 's/$/\tSynechocystis/'
wget -qO - http://rest.kegg.jp/link/ec/hpse | sed -e 's/$/\tHydrogenophaga/'
) > data/KEGGgene_EC_organism.tab

# Module descriptions
wget -qO data/module_description.tab http://rest.kegg.jp/list/module

# Module to compound connection
wget -qO data/module_compound.tab http://rest.kegg.jp/link/compound/module

# Module to compound full
wget -qO data/module_reaction.tab http://rest.kegg.jp/link/reaction/module
wget -qO data/reaction_compound.tab http://rest.kegg.jp/link/compound/reaction

# EC descriptions
wget -qO data/EC_description.tab http://rest.kegg.jp/list/ec
