library(argparse)
library(monocle)
library(stringr)
library(dplyr)
library(Matrix)
library(tidyr)


source("custom_monocle_da_functions.R")


get_gene_vs_gRNAgroup_1Mb_DE <- function(gRNA_group, cds_obj, gene_gRNAgroup_pair_table){
  
  genes_to_test <- subset(gene_gRNAgroup_pair_table, gRNAgroup == gRNA_group)
  
  #print(head(genes_to_test))
  
  ENSG_list <- as.character(genes_to_test$ENSG.targetgene) #genes to test against gRNA group, called within 1 Mb
  
  gene_cds <- cds_obj[rownames(cds_obj) %in% ENSG_list, ] #wont require all the ensgs are in the list
  #print("Size_Factor" %in% colnames(pData(gene_cds)))
  #print(head(pData(gene_cds)))
  
  pData(gene_cds)$gRNA_detected <- grepl(paste0('(^|_)', gRNA_group, '(_|$)'), pData(gene_cds)$gene) #this should refer to the true/false input if a gRNA is detected, esp as i re-run with gRNAgroup
  print(table(pData(gene_cds)$gRNA_detected))
  print(gRNA_group)
  gene_cds <- gene_cds[, !is.na(pData(gene_cds)$gRNA_detected)]
  #print("Size_Factor" %in% colnames(pData(gene_cds)))
  #print("No gRNA detected: ")
  #print( w_NA -  no_NA) #print the number of cells w/o gRNA detected
  
  #print(table(pData(gene_cds)$gRNA_detected))
  
  #print("About to start DE")
  #print(dim(gene_cds)) 
  #closeAllConnections()

  DE_results <- myDifferentialGeneTest(gene_cds, fullModelFormulaStr = "~gRNA_detected", cores=1)
  print(head(DE_results)) 
  #print(DE_results[1,'qval'])
  #tried to merge but i think that's confusing
  #merged_df <- merge(gene_gRNAgroup_pair_table, DE_results, by.x = "ENSG.targetgene", by.y = "id")
  
  #DE_results$gRNA_group <- gRNA_group
  
  #print(head(DE_results))
  #print(head(gene_gRNAgroup_pair_table))
  DE_results$gRNA_group <- gRNA_group
  
  #ENSG_list HAS THE GENE SOMEHOW
  

  ####EXPORT THE MEANS
  
  #print(table(pData(gene_cds)$gRNA_detected))
  
  yes_gRNA_cds <- gene_cds[,pData(gene_cds)$gRNA_detected]
  no_gRNA_cds <- gene_cds[,!pData(gene_cds)$gRNA_detected]
  #print("Size_Factor" %in% colnames(pData(yes_gRNA_cds)))
  #print("Size_Factor" %in% colnames(pData(no_gRNA_cds)))
  print(yes_gRNA_cds)
  print(no_gRNA_cds)
  #print('1')
  yes_gRNA_rowmean <- Matrix::rowSums(t( t(exprs(yes_gRNA_cds)) / pData(yes_gRNA_cds)$Size_Factor) ) / ncol(yes_gRNA_cds) #have to double transpose to make sure the sizefactor will divide the right way, transpose again to get it to work
  #print('2')
  #print(no_gRNA_cds)
  no_gRNA_rowmean <- Matrix::rowSums(t( t(exprs(no_gRNA_cds)) / pData(no_gRNA_cds)$Size_Factor) ) / ncol(no_gRNA_cds)

  DE_results$yes_gRNA_rowmean <- yes_gRNA_rowmean
  DE_results$no_gRNA_rowmean <- no_gRNA_rowmean

  return(DE_results)
  
  
}

###TROUBLESHOOTING###


#source("~/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/Vol10_PROJECTS/170803_CRISPRQTL/bin/custom_monocle_da_functions.R")


#test.rds <- readRDS("~/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/Vol10_PROJECTS/170803_CRISPRQTL/results/171012_K1000_analysis/K1000.NOCHIM.sizefact.dispersions.exprs.rds")

#test.gene_gRNAgroup_pair_table <- read.delim("~/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/Vol10_PROJECTS/170803_CRISPRQTL/results/171012_K1000_analysis/171018_gRNAxGene_association_tests/171018_gRAx2Mbgene_PAIRS.txt")

#test.gRNA_group <- "ACTB_TSS"

#undebug(get_gene_vs_gRNAgroup_1Mb_DE)

#test.result <- get_gene_vs_gRNAgroup_1Mb_DE(test.gRNA_group, test.rds, test.gene_gRNAgroup_pair_table)

###TROUBLESHOOTING END ###




parser=argparse::ArgumentParser(description = "Script to perform DE on gRNAgroups")
parser$add_argument("cds_obj")
parser$add_argument("gRNA_group")
parser$add_argument("gene_gRNAgroup_pair_table")
parser$add_argument("output")

args = parser$parse_args()

K1000_cds <- readRDS(args$cds_obj) #43000 cells

READIN.gene_gRNAgroup_pair_table <- as.data.frame(readr::read_delim(args$gene_gRNAgroup_pair_table, delim='\t'))

DE_results <- get_gene_vs_gRNAgroup_1Mb_DE(args$gRNA_group, K1000_cds, READIN.gene_gRNAgroup_pair_table)

write.table(DE_results, file = args$output , quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
