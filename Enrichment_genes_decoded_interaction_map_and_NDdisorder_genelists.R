# ENRICHMENT SCRIPT 

# this code tests significant overlap between the list previously obtained (genes decoded) and
# genes associated to neurodevelopmental disorders.


library(ggplot2)
library(readxl)

data_path='NDdisorders_genelist/'

# LOADING OF GENE DECODING DATA

# source main enrichment script, edit /path/genelistOverlap.r
source("genelistOverlap.r")

# rootpath where files are located, edit the path of the root folder
rootpath = "00_Genes_decoded_interaction_map_sexandvolsretained_regressed_p005/"

# this reads the gene list, edit the folder name (path) and the list name (txt file) 
postimg_fdr05 = unique(as.character(read.delim(file.path(rootpath,"genesample_pos_thresh_results.txt"),header=FALSE)$V1))

#this reads the list of genes highly expressed in the brain 
high_brain_expr = unique(as.character(read.delim(file.path("BrainExpressed_Gandal/Gandal_Nature2022_16796_genes_highly_expressed_in_the_brain.txt"), header = FALSE)$V1))


# background totals to use for enrichment analyses 
high_expr_brain_background = as.numeric(length(high_brain_expr))

### Now do a pruning of the decoded gens - as to include only genes that are highly expressed in the brain
postimg_fdr05_filtered = intersect(postimg_fdr05, high_brain_expr)


output_dir='03_ENRICHMENT_between_genes_decoded_interaction_map_and_NDdisorder_genelists/'
# this creates the output directory
if (!file.exists(output_dir)) {
  dir.create(output_dir)
}



#---------------------------------------------------------------------------------------------------------------------------------------------
## ASD-related enrichment analyses with decoded genes
# pre-allocate array with gene lists to use
gene_lists2use = c("SFARI",
                   "Rare inherited ASD variants (Wilfert, Nature Genetics, 2021)",
                   "ASD (Satterstrom, Cell, 2020)",
                   "ASD up or down regulated genes, whole cortex (Gandal, Nature, 2022)",
                   "ASD upregulated genes, whole cortex (Gandal, Nature, 2022)",
                   "ASD downregulated genes, whole cortex (Gandal, Nature, 2022)",
                   "ASD common variants (Grove, Nat. Genet., 2019)",
                   "Bipolar Disorder (Gandal, Science, 2018)",
                   "Schizophrenia (Gandal, Science, 2018)",
                   "ADHD genes (Satterstrom et al., Nat. Neurosci., 2019)"
)


gene_lists_savename = c("SFARI",
                        "Rare_inherited_ASD variants_Wilfert",
                        "ASD_Satterstrom_2020",
                        "ASD_up_or_down_regulated_genes_whole_cortex_Gandal",
                        "ASD_up_regulated_genes_whole_cortex_Gandal",
                        "ASD_down_regulated_genes_whole_cortex_Gandal",
                        "ASD_common_variants_Grove",
                        "Bipolar_Disorder",
                        "Schizophrenia",
                        "ADHD_genes_Satterstrom"
)

cols2use = c("genelist","OR","pval","fdr")
gsea_res = data.frame(matrix(nrow = length(gene_lists2use), ncol = length(cols2use)))
rownames(gsea_res) = gene_lists2use
colnames(gsea_res) = cols2use



# load SFARI genes
sfari = read.csv(file.path(data_path,"SFARI-Gene_genes_04-03-2023release_04-19-2023export.csv"))
sfari_genes = sfari$gene_symbol
sfari_genes = sfari_genes[is.element(sfari_genes, high_brain_expr)]

# load Wilfert et al., 2021 genes
wilfert = read.csv(file.path(data_path,"wilfert2021_supptab17_privateInheritedLGD.csv"))
wilfert_genes = wilfert$Gene
wilfert_genes = wilfert_genes[is.element(wilfert_genes, high_brain_expr)]

# load Satterstrom asd genes
satterstrom = read.csv(file.path(data_path, "satterstrom_suppT2.csv"))
satterstrom_genes = satterstrom$gene
satterstrom_genes = satterstrom_genes[is.element(satterstrom_genes, high_brain_expr)]

# load Gandal all genes 2022
gandal = read.csv(file.path(data_path,"gandal_2022_nature_suppdata03.csv"))
mask = gandal$WholeCortex_ASD_FDR<=0.05
gandal_all_genes = gandal$external_gene_name[mask]
mask = gandal_all_genes==""
gandal_all_genes = gandal_all_genes[!mask]
gandal_all_genes = gandal_all_genes[is.element(gandal_all_genes, high_brain_expr)]

# load Gandal whole cortex upregulated
mask = (gandal$WholeCortex_ASD_FDR<=0.05) & (gandal$WholeCortex_ASD_logFC>0)
gandal_all_up_genes = gandal$external_gene_name[mask]
mask = gandal_all_up_genes==""
gandal_all_up_genes = gandal_all_up_genes[!mask]
gandal_all_up_genes = gandal_all_up_genes[is.element(gandal_all_up_genes, high_brain_expr)]

# load Gandal whole cortex downregulated
mask = (gandal$WholeCortex_ASD_FDR<=0.05) & (gandal$WholeCortex_ASD_logFC<0)
gandal_all_down_genes = gandal$external_gene_name[mask]
mask = gandal_all_down_genes==""
gandal_all_down_genes = gandal_all_down_genes[!mask]
gandal_all_down_genes = gandal_all_down_genes[is.element(gandal_all_down_genes, high_brain_expr)]

# load Grove asd genes
grove = read.csv(file.path(data_path, "genelist_autism_Grove_NatureGenetics2019_Supplementary_Table4.csv"))
grove_genes = grove$Gene
grove_genes = grove_genes[is.element(grove_genes, high_brain_expr)]

# load Gandal - Bipolar 
gandal2018 = read.csv(file.path(data_path,"gandal2018_degenes.csv"))
mask = gandal2018$BD.fdr<=0.05
bd_genes = gandal2018$gene_name[mask]
mask = bd_genes==""
bd_genes = bd_genes[!mask]
bd_genes = bd_genes[is.element(bd_genes, high_brain_expr)]

# load Gandal - Schizophrenia 
gandal2018 = read.csv(file.path(data_path,"gandal2018_degenes.csv"))
mask = gandal2018$SCZ.fdr<=0.05
scz_genes = gandal2018$gene_name[mask]
mask = scz_genes==""
scz_genes = scz_genes[!mask]
scz_genes = scz_genes[is.element(scz_genes, high_brain_expr)]

# load Satterstrom - ADHD  
adhd_satterstrom = read.csv(file.path(data_path,"genelist_Sattrestrom_gene_mutated_in_ADHD.csv"))
adhd_satterstrom = adhd_satterstrom$Gene
adhd_satterstrom_genes = adhd_satterstrom[is.element(adhd_satterstrom, high_brain_expr)]


######### gene enrichment ######### 

gene_lists4loop = list(sfari_genes,
                       wilfert_genes,
                       satterstrom_genes,
                       gandal_all_genes,
                       gandal_all_up_genes,
                       gandal_all_down_genes,
                       grove_genes,
                       bd_genes,
                       scz_genes,
                       adhd_satterstrom_genes
                       )



list_name2use = "group_table"

result_path=output_dir


# Enrichment of all lists
for (ilist in 1:length(gene_lists4loop)){
  list1 = as.data.frame(postimg_fdr05_filtered)
  list2 = gene_lists4loop[[ilist]]
  backgroundTotal = high_expr_brain_background
  overlap_res = genelistOverlap(list1 = list1, list2 = as.data.frame(list2),
                                backgroundTotal = backgroundTotal, print_result = TRUE) 
  
 # write.csv(overlap_res[[1]]$overlapping_genes, file = file.path(paste0(result_path, gene_lists_savename[ilist],".csv")))
  
 
  gsea_res[gene_lists2use[ilist],"genelist"] = gene_lists2use[ilist]
  gsea_res[gene_lists2use[ilist],"genelist length"] = length(list2)
  gsea_res[gene_lists2use[ilist],"OR"] = overlap_res[[1]]$OR
  gsea_res[gene_lists2use[ilist],"pval"] = overlap_res[[1]]$hypergeo_p
  gsea_res[gene_lists2use[ilist],"Gene Overlap"] = length(overlap_res[[1]]$overlapping_genes)
  
  
  
}
gsea_res[,"fdr"] = p.adjust(gsea_res[,"pval"], method = "fdr")
gsea_res
write.csv(gsea_res, file=paste0(output_dir,'genes_decoded_interaction_map_and_NDdisorders_genelists.csv'))

