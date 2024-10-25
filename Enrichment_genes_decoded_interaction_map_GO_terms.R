# ENRICHMENT SCRIPT 

# this code tests significant overlap between the list previously obtained (genes decoded) and
# GO Terms. 

library(ggplot2)
library(readxl)


# source main enrichment script, edit /path/genelistOverlap.r
source("genelistOverlap.r")

# rootpath where files are located, edit the path of the root folder
rootpath = "../00_Genes_decoded_interaction_map_sexandvolsretained_regressed_p005"


output_dir='../01_ENRICHMENT_between_genes_decoded_interaction_map_GOTerms/'
# this creates the output directory
if (!file.exists(output_dir)) {
  dir.create(output_dir)
}

######## BACKGROUND GENELIST ######## 

#this reads the list of genes highly expressed in the brain 
high_brain_expr = unique(as.character(read.delim(file.path("../BrainExpressed_Gandal/Gandal_Nature2022_16796_genes_highly_expressed_in_the_brain.txt"), header = FALSE)$V1))

# background totals to use for enrichment analyses 
high_expr_brain_background = as.numeric(length(high_brain_expr))                   

# this reads in the synaptic genes
synaptic_genes_full_list <-  unique(as.character(read.delim(file.path("../Synaptic_genes/syngo_genes_fixed.txt"), header = FALSE)$V1))
          
     
######### GENELIST 1: the list of genes in positively correlated with interaction maps ##############

# this reads the gene list, edit the folder name (path) and the list name (txt file) 
postimg_fdr05 = unique(as.character(read.delim(file.path(rootpath,"genesample_pos_thresh_results.txt"),header=FALSE)$V1))
#  Now do a pruning of the decoded gens - as to include only genes that are highly expressed in the brain
postimg_fdr05_filtered = intersect(postimg_fdr05, high_brain_expr)



######### List of genes to run the enrichment with:

GO_lists = list.files(path="../Additional_checks_GOTerms/")
index=1
all_lists=list()

for (i in GO_lists){
 # print(i)
  filename=file.path("../Additional_checks_GOTerms", i)
  print(filename)
  all_lists[[index]] <-  unique(as.character(read.delim(filename, header = FALSE)$V1))
  
  index=index+1
  
}

# saving all GO terms lists in a human readable way, so that they can be saved and shared
# List of column names
column_names <- gene_lists2use

# Find the maximum length among the list elements
max_length <- max(sapply(all_lists, length))

# Convert to data frame, filling shorter vectors with NAs
my_df <- as.data.frame(lapply(all_lists, function(x) {
  length(x) <- max_length  # Extend shorter vectors
  x[is.na(x)] <- ""        # Replace NAs with blank spaces
  
  return(x)
}))

# Set the column names
colnames(my_df) <- column_names

# View the data frame
write.csv(my_df, "GOTerms_lists_Human_readable_combined.csv")


# pre-allocate array with gene lists to use
gene_lists2use = GO_lists


gene_lists2use <- gsub("_", " ",   gene_lists2use)
gene_lists2use <- gsub(".txt", "",   gene_lists2use)
gene_lists2use <-sub("(.{9})(.*)", "\\1:\\2", gene_lists2use)

cols2use = c("genelist", "genelist length","OR","pval","fdr","Gene Overlap")
gsea_res = data.frame(matrix(nrow = length(gene_lists2use), ncol = length(cols2use)))
rownames(gsea_res) = gene_lists2use
colnames(gsea_res) = cols2use

gene_lists4loop = all_lists

list_name2use = "group_table"


result_path=output_dir


# Enrichment of all lists 
for (ilist in 1:length(gene_lists4loop)){
  list1 = as.data.frame(postimg_fdr05_filtered)
  list2 = gene_lists4loop[[ilist]]
  list2 = as.data.frame(list2[is.element(list2,high_brain_expr)])
  
  backgroundTotal = high_expr_brain_background
  overlap_res = genelistOverlap(list1 = list1, list2 = as.data.frame(list2),
                                backgroundTotal = backgroundTotal, print_result = TRUE) 
  
  
  gsea_res[gene_lists2use[ilist],"genelist"] = gene_lists2use[ilist]
  gsea_res[gene_lists2use[ilist],"genelist length"] =  dim(list2)[1]
  gsea_res[gene_lists2use[ilist],"OR"] = overlap_res[[1]]$OR
  gsea_res[gene_lists2use[ilist],"pval"] = overlap_res[[1]]$hypergeo_p
  gsea_res[gene_lists2use[ilist],"Gene Overlap"] = length(overlap_res[[1]]$overlapping_genes)
  
  # write.table(overlap_res[[1]]$overlapping_genes, file = paste0(result_path, gene_lists2use[ilist], "OVERLAP with genes decoded.txt"), 
          #     sep = "\n", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
}
gsea_res[,"fdr"] = p.adjust(gsea_res[,"pval"], method = "fdr")
gsea_res
write.csv(gsea_res, file = file.path(paste0(result_path,"output_enrichment.csv")))





