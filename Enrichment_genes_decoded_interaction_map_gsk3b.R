# ENRICHMENT SCRIPT 

# this code tests significant overlap between the list previously obtained (genes decoded) and
# the GSK3B interactome. 

library(ggplot2)
library(readxl)


# source main enrichment script, edit /path/genelistOverlap.r
source("genelistOverlap.r")

# rootpath where files are located, edit the path of the root folder
rootpath = "00_Genes_decoded_interaction_map_sexandvolsretained_regressed_p005/"



output_dir='02_ENRICHMENT_between_genes_decoded_interaction_map_and_GSK3B/'
# this creates the output directory
if (!file.exists(output_dir)) {
  dir.create(output_dir)
}

######## BACKGROUND GENELIST ######## 

#this reads the list of genes highly expressed in the brain 
high_brain_expr = unique(as.character(read.delim(file.path("BrainExpressed_Gandal/Gandal_Nature2022_16796_genes_highly_expressed_in_the_brain.txt"), header = FALSE)$V1))


# background totals to use for enrichment analyses 
high_expr_brain_background = as.numeric(length(high_brain_expr))                  


######### GENELIST 1: the list of genes in positively correlated with interaction maps ##############

# this reads the gene list. edit the folder name (path) and the list name (txt file) 
postimg_fdr05 = unique(as.character(read.delim(file.path(rootpath,"genesample_pos_thresh_results.txt"),header=FALSE)$V1))
#  Now do a pruning of the decoded gens - as to include only genes that are highly expressed in the brain
postimg_fdr05_filtered = intersect(postimg_fdr05, high_brain_expr)



######### list of genes to run the enrichment with:
gsk3b_medium = unique(as.character(read.delim("Interactomes/gsk3b_medium_confidence_maxinteraction500.txt",header=FALSE)$V1))
gsk3b_medium_filtered = intersect(gsk3b_medium, high_brain_expr)

gsk3b_medium_synaptic <- unique(as.character(read.delim("Interactomes/gsk3b_synaptic_interactome_as_downloaded_from_Syngo.txt",header=FALSE)$V1))
gsk3b_medium_synaptic_brainexpr <- intersect(gsk3b_medium_synaptic,high_brain_expr)


# pre-allocate array with gene lists to use
gene_lists2use = c("GSK3B interactors",
                   "Synaptic GSKB interactors"
)

cols2use = c("genelist", "genelist length","OR","pval","fdr","Gene Overlap")
gsea_res = data.frame(matrix(nrow = length(gene_lists2use), ncol = length(cols2use)))
rownames(gsea_res) = gene_lists2use
colnames(gsea_res) = cols2use

gene_lists4loop = list(gsk3b_medium_filtered,
                       gsk3b_medium_synaptic_brainexpr
)


list_name2use = "group_table"

result_path=output_dir


# Enrichment GSK3b 
for (ilist in 1:length(gene_lists4loop)){
  list1 = as.data.frame(postimg_fdr05_filtered)
  list2 = gene_lists4loop[[ilist]]
  backgroundTotal = high_expr_brain_background
  overlap_res = genelistOverlap(list1 = list1, list2 = as.data.frame(list2),
                                backgroundTotal = backgroundTotal, print_result = TRUE) 
  
  
  gsea_res[gene_lists2use[ilist],"genelist"] = gene_lists2use[ilist]
  gsea_res[gene_lists2use[ilist],"genelist length"] = length(list2)
  gsea_res[gene_lists2use[ilist],"OR"] = overlap_res[[1]]$OR
  gsea_res[gene_lists2use[ilist],"pval"] = overlap_res[[1]]$hypergeo_p
  gsea_res[gene_lists2use[ilist],"Gene Overlap"] = length(overlap_res[[1]]$overlapping_genes)
  
   write.table(overlap_res[[1]]$overlapping_genes, file = paste0(result_path, gene_lists2use[ilist], "OVERLAP with genes decoded.txt"), 
               sep = "\n", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
}
gsea_res[,"fdr"] = p.adjust(gsea_res[,"pval"], method = "fdr")
gsea_res
write.csv(gsea_res, file = file.path(paste0(result_path,"output_enrichment.csv")))




