library(argparse)
library(monocle)
library(stringr)
library(dplyr)
library(Matrix)
library(tidyr)

source("custom_monocle_da_functions.R")


get_gene_vs_gRNAgroup_1Mb_DE <- function(gRNA_group, cds_obj, gene_gRNAgroup_pair_table, reference_cells_rds=NULL, guide=NULL){
  
  genes_to_test <- subset(gene_gRNAgroup_pair_table, gRNAgroup == gRNA_group)
  
  #print(head(genes_to_test))
  
  ENSG_list <- as.character(genes_to_test$ENSG.targetgene) #genes to test against gRNA group, called within 1 Mb
  
  gene_cds <- cds_obj[rownames(cds_obj) %in% ENSG_list, ] #wont require all the ensgs are in the list
  #print("Size_Factor" %in% colnames(pData(gene_cds)))
  #print(head(pData(gene_cds)))
  # Can test either the whole group or a specific guide
  if (is.null(guide)) {
    pData(gene_cds)$gRNA_detected <- grepl(paste0('(^|_)', gRNA_group, '(_|$)'), pData(gene_cds)$gene) #this should refer to the true/false input if a gRNA is detected, esp as i re-run with gRNAgroup
  } else {
    pData(gene_cds)$gRNA_detected <- grepl(paste0('(^|_)', guide, '(_|$)'), pData(gene_cds)$barcode
)
  }

  if (!is.null(reference_cells_rds)) {
    reference_cells = colnames(gene_cds) %in% readRDS(reference_cells_rds)
    positive_cells = !is.na(pData(gene_cds)$gRNA_detected) & pData(gene_cds)$gRNA_detected 
    gene_cds = gene_cds[, reference_cells | positive_cells]
  }

  gene_cds <- gene_cds[, !is.na(pData(gene_cds)$gRNA_detected)]
  gene_cds$batch = with(pData(gene_cds), paste(prep_batch, within_batch_chip, within_chip_lane, sep='.'))
  DE_results <- myDifferentialGeneTest(gene_cds, fullModelFormulaStr = "~gRNA_detected+percent.mito+guide_count+prep_batch", reducedModelFormulaStr='~percent.mito+guide_count+prep_batch', cores=1, verbose=TRUE)
  DE_results$gRNA_group <- gRNA_group
 
  if (!is.null(guide)) {
    DE_results$guide = guide
  } 
  
  yes_gRNA_cds <- gene_cds[,pData(gene_cds)$gRNA_detected]
  no_gRNA_cds <- gene_cds[,!pData(gene_cds)$gRNA_detected]
  yes_gRNA_rowmean <- Matrix::rowMeans(t( t(exprs(yes_gRNA_cds)) / pData(yes_gRNA_cds)$Size_Factor))[DE_results$id] #have to double transpose to make sure the sizefactor will divide the right way, transpose again to get it to work
  no_gRNA_rowmean <- Matrix::rowMeans(t( t(exprs(no_gRNA_cds)) / pData(no_gRNA_cds)$Size_Factor))[DE_results$id]

  DE_results$yes_gRNA_rowmean <- yes_gRNA_rowmean
  DE_results$no_gRNA_rowmean <- no_gRNA_rowmean

  return(DE_results)
}

parser=argparse::ArgumentParser(description = "Script to perform DE on gRNAgroups")
parser$add_argument("cds_obj")
parser$add_argument("gRNA_group")
parser$add_argument("gene_gRNAgroup_pair_table")
parser$add_argument("output")
parser$add_argument('--reference_cells_rds', default=NULL)
parser$add_argument('--guide', default=NULL, help='Optional can also specify a specific guide to test individually.')

args = parser$parse_args()

K1000_cds <- readRDS(args$cds_obj) #43000 cells

READIN.gene_gRNAgroup_pair_table <- as.data.frame(readr::read_delim(args$gene_gRNAgroup_pair_table, delim='\t'))

DE_results <- get_gene_vs_gRNAgroup_1Mb_DE(args$gRNA_group, K1000_cds, READIN.gene_gRNAgroup_pair_table, args$reference_cells_rds, guide=args$guide)

write.table(DE_results, file = args$output , quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
