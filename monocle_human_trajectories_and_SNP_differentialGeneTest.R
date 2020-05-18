
library(reticulate)
library("devtools")
library(DESeq)
library(plyr)
library(dplyr)
library(stringr)
library("VariantAnnotation")
suppressMessages(library(tidyr))
suppressMessages(library(Biobase))
library(knitr)
library(reshape2)
library(ggplot2)
library(Geno1000matrixSingleCell)
library(monocle)
library(radmixture)
library(randomcoloR)

g_1000<-read.delim("/Users/weiz/Downloads/recoded_1000G.raw",header = T,sep = " ")

g_1000<-read.delim("/Users/weiz/Downloads/recoded_1000G.raw",header = T,sep = " ")
rownames(g_1000) <- g_1000[,2]
g_1000[1:10,1:10]
human_list<-data.frame(rownames(g_1000))

write.table(human_list,"/Users/weiz/Desktop/R_code/e11/human_list.txt",quote = F,sep = "\t")

unprocessed_single_cell_data_total_transposed <-g_1000[,7:39188]
population<- read.delim("/Users/weiz/Downloads/igsr_samples.tsv",header = T)



datadata<-as.data.frame(as.character(colnames(unprocessed_single_cell_data_total_transposed)) )
mutation_genotype<-datadata %>% separate("as.character(colnames(unprocessed_single_cell_data_total_transposed))", into=c("ID", "genotype"), "_")
#mutation_genotype <-mutation_genotype[7:length(mutation_genotype[,1]),]

mutation_genotype<-cbind(mutation_genotype,mutation_genotype$ID,mutation_genotype$genotype)
write.table(mutation_genotype,"/Users/weiz/Desktop/R_code/e11/mutation_genotype.txt",quote = F,sep = "\t")
#area<-"EAS"
pdf("plot_monocle_human_areaEAS.pdf",20,20)
human_monocle<-function(area){
  
  
  population<-as.data.frame(population[which(population$Superpopulation.code==area),])
  
  area_unprocessed_single_cell_data_total_transposed<-intersect(population$Sample.name,rownames(unprocessed_single_cell_data_total_transposed))
  
  unprocessed_single_cell_data_total_transposed_area<-unprocessed_single_cell_data_total_transposed[area_unprocessed_single_cell_data_total_transposed,]
  
  
  
  Geno1000matrix_expr_matrix<-t(unprocessed_single_cell_data_total_transposed_area)
  Geno1000matrix_expr_matrix[is.na(Geno1000matrix_expr_matrix)] <- 0
  rpc_matrix <- relative2abs(Geno1000matrix, method = "tpm_fraction")
  Geno1000matrix <- newCellDataSet(as.matrix(rpc_matrix))
  Geno1000matrix <- estimateSizeFactors(Geno1000matrix)
  Geno1000matrix <- estimateDispersions(Geno1000matrix)
  Geno1000matrix <- detectGenes(Geno1000matrix, min_expr = 0.1)
  print(head(fData(Geno1000matrix)))
  expressed_genes <- row.names(subset(fData(Geno1000matrix), num_cells_expressed >= 10))
  length(expressed_genes)
  print(head(pData(Geno1000matrix)))
  
  Geno1000matrix <- detectGenes(Geno1000matrix, min_expr = 0.1)
  
  L <- log(exprs(Geno1000matrix[expressed_genes,]))
  
  melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))
  
  qplot(value, geom="density", data=melted_dens_df) +  stat_function(fun = dnorm, size=0.5, color='red') +
    xlab("Standardized log(FPKM)") +
    ylab("Density")
  
  
  disp_table <- dispersionTable(Geno1000matrix)
  head(disp_table)
  
  
  unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
  Geno1000matrix <- setOrderingFilter(Geno1000matrix, unsup_clustering_genes$gene_id)
  print(plot_ordering_genes(Geno1000matrix))
  plot_pc_variance_explained(Geno1000matrix, return_all = F) # norm_method = 'log',
  Geno1000matrix <- reduceDimension(Geno1000matrix, max_components=2, num_dim = 6,
                          reduction_method = 'tSNE', verbose = T)
  Geno1000matrix <- clusterCells(Geno1000matrix, num_clusters=2)
  head(fData(Geno1000matrix))
  print(plot_cell_clusters(matrix_resample, 1, 2,color_by=population_names))
  matrix_resample <- estimateDispersions(Geno1000matrix,method='blind',sharingMode="fit-only")
  
  disp_table <- dispersionTable(matrix_resample)
  ordering_genes <- subset(disp_table,
                           mean_expression >= 0.01 &
                             dispersion_empirical >= 1 * dispersion_fit)$gene_id
  
  
  matrix_resample <- setOrderingFilter(matrix_resample, ordering_genes)
  print(plot_ordering_genes(matrix_resample))
  
  
  #matrix_resample <- reduceDimension(matrix_resample, max_components=2,reduction_method = 'DDRTree', verbose = T)
  matrix_resample <- orderCells(matrix_resample)
  print(plot_cell_trajectory(matrix_resample, color_by="State"))
  
  state_file<-plot_cell_trajectory(matrix_resample, color_by="State")
  print(plot_cell_trajectory(matrix_resample, color_by = "Pseudotime"))
  
  Pseudotime_file<-plot_cell_trajectory(matrix_resample, color_by = "Pseudotime")
  
  
  Pseudotime_data_table<-cbind (sample_name=state_file[["data"]][["sample_name"]],dim1=state_file[["data"]][["data_dim_1"]],dim2=state_file[["data"]][["data_dim_2"]],State=state_file[["data"]][["State"]],Pseudotime=Pseudotime_file[["data"]][["Pseudotime"]])
  
  write.table( Pseudotime_data_table,paste("Pseudotime_data_table",area,sep = "_",".txt"),quote = F,sep="\t")
  
  
  
  position_aeas <- match(colnames(matrix_resample),population$Sample.name)
  
  
  population_names<-population[position_aeas,5]
  
  
  pdf("plot_monocle_human_area_State.pdf")
  
  print(plot_cell_trajectory(matrix_resample, color_by="State"))
  
  dev.off()
  pdf("plot_monocle_human_area_Pseudotime.pdf")
  print(plot_cell_trajectory(matrix_resample, color_by = "Pseudotime"))
  dev.off()
  pdf("plot_monocle_human_area_population_names.pdf")
  print(plot_cell_trajectory(matrix_resample, color_by=population_names) )
  
  dev.off()
  
  pdf("group_plot_monocle_human_area_population_names.pdf")
  print(plot_cell_clusters(matrix_resample, 1, 2,color_by=population_names))
  dev.off()
}
dev.off()

unique(population_names)

pdf("plot_monocle_human_area3.pdf",20,20)

for (area in unique(population$Superpopulation.code))
  
{
  human_monocle(area)
  
}

dev.off()


##############
matrix_resample <- estimateDispersions(Geno1000matrix)

disp_table <- dispersionTable(matrix_resample)
ordering_genes <- subset(disp_table,
                         mean_expression >= 0.1 &
                           dispersion_empirical >= 0.0001 * dispersion_fit)$gene_id


matrix_resample <- setOrderingFilter(matrix_resample, ordering_genes)
plot_ordering_genes(matrix_resample)

matrix_resample <- clusterCells(matrix_resample, verbose = F)



diff_test_res <- differentialGeneTest(matrix_resample[expressed_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, pval < 0.01))

datadata<-as.data.frame(as.character(sig_gene_names) )
mutation_genotype<-datadata %>% separate("as.character(sig_gene_names)", into=c("ID", "genotype"), "_")
write.table(mutation_genotype,"mutation_genotype_001.txt",quote = F,sep = "\t")



library(biomaRt)


mart.snp <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "grch37.ensembl.org")

dbsnp = useMart("snp", dataset = "hsapiens_snp")
listMarts
getENSG <- function(rs = "rs3043732", mart = mart.snp) {
  results <- getBM(attributes = c("refsnp_id", "ensembl_gene_stable_id"),
                   filters    = "snp_filter", values = rs, mart = mart)
  return(results)
}

# default parameters
getENSG()
refsnp_id ensembl_gene_stable_id
1 rs3043732        ENSG00000175445

# or supply rs ID
getENSG(rs = "rs224550")




bks <- seq(-3.1, 3.1, by = 0.1)
hmcols <- blue2green2red(length(bks) - 1)
pdf(paste("Human trajectory plot ofplot_siggene_each_plot_pseudotime_heatmap001.pdf",sep=""))

print(plot_pseudotime_heatmap(matrix_resample[sig_gene_names,],
                              num_clusters = 8,
                              cores = 1,
                              show_rownames = T))

dev.off()


expressed_genes <- row.names(subset(fData(matrix_resample), num_cells_expressed >= 400))
diff_test_res <- differentialGeneTest(matrix_resample[expressed_genes,],fullModelFormulaStr = "~clusters")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
matrix_resample <- setOrderingFilter(matrix_resample, ordering_genes)
plot_ordering_genes(matrix_resample)
matrix_resample <- reduceDimension(matrix_resample, max_components = 2)
matrix_resample <- orderCells(matrix_resample)

