# Copy over scripts from original directories
cp /net/trapnell/vol1/ajh24/proj/2017crispr_qtl/data/reads/2017_11_10_K1000.high_moi/171018_gRNAx2Mbgene_DE_means.R get_deg.pilot.R
 
cp /net/trapnell/vol1/ajh24/proj/2017crispr_qtl/data/reads/2018_04_24_K6000/171018_gRNAx2Mbgene_DE_means.R get_deg.at_scale.R

# Copy over the custom monocle code
cp ~gasperim/PROJECTS/170803_CRISPRQTL/bin/custom_monocle_da_functions.R .

#Copy over the 50k reference cell set used for the at-scale screen DEG calling
cp /net/trapnell/vol1/ajh24/proj/2017crispr_qtl/data/reads/2018_04_24_K6000/reference_cells.50K.rds .

mv reference_cells.50K.rds 50k_reference_cells.rds

