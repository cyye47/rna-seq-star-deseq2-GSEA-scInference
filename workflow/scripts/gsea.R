# --- Load Required Packages ---
suppressPackageStartupMessages({
  #library(clusterProfiler)
  library(org.Sc.sgd.db)
  #library(enrichplot)
  library(DOSE)
  library(ggplot2)
  library(readr)
  library(dplyr)
})

# # --- Input and Output Directories ---
# args <- commandArgs(trailingOnly = TRUE)
# input_file <- args[1]
# print(input_file)
# out_dir <- args[2]
# print(out_dir)

out_dir <- snakemake@output[["output_dir"]]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Read DESeq2 Results ---
df <- read_tsv(snakemake@input[["diffexpr"]])

# Expecting columns: gene, log2FoldChange
df <- df %>% filter(!is.na(log2FoldChange))

# --- Convert Gene Symbols (ORFs) to Entrez IDs ---
gene_df <- bitr(df$gene, fromType = "ORF",
                toType = "ENTREZID",
                OrgDb = org.Sc.sgd.db)

# --- Merge Data and Rank Genes ---
merged_df <- inner_join(df, gene_df, by = c("gene" = "ORF"))
gene_list <- merged_df$log2FoldChange
names(gene_list) <- merged_df$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

# --- Functions for pathway analysis ---
# --- (GO) ---
gsea_go <- function(gene_list) {
    gseGO(
        geneList = gene_list,
        OrgDb = org.Sc.sgd.db,
        ont = "BP",
        nPerm = 1000,
        minGSSize = 10,
        pvalueCutoff = 0.05,
        verbose = FALSE
        )
}

# --- (KEGG) ---
gsea_kegg <- function(gene_list) {
    gseKEGG(
        geneList = gene_list,
        organism = 'sce',  # KEGG code for S. cerevisiae
        nPerm = 1000,
        minGSSize = 10,
        pvalueCutoff = 0.05,
        verbose = FALSE
        )

}


##############################
# For all disregulated genes #
##############################

gsea_go_result <- gsea_go(gene_list)
gsea_kegg_result <- gsea_kegg(gene_list)

gsea_go_file <- file.path(out_dir, "gsea_go_dge_all_results.tsv")
gsea_kegg_file <- file.path(out_dir, "gsea_kegg_dge_all_results.tsv")

if(nrow(gsea_go_result@result) > 0) {

  write_tsv(gsea_go_result@result, gsea_go_file)

  ## GO Dotplot
  pdf(file.path(out_dir, "gsea_go_dge_all_dotplot.pdf"), width = 8, height = 6)
  dotplot(gsea_go_result, showCategory = 20)
  dev.off()
  
  ## GO Ridgeplot
  pdf(file.path(out_dir, "gsea_go_dge_all_ridgeplot.pdf"), width = 8, height = 6)
  ridgeplot(gsea_go_result, showCategory = 20)
  dev.off()

  ## GO Top Pathway
  top_go_id <- gsea_go_result@result$ID[1]
  pdf(file.path(out_dir, "gsea_go_dge_all_top_pathway.pdf"), width = 8, height = 6)
  gseaplot2(gsea_go_result, geneSetID = top_go_id)
  dev.off()

} else {

  write_tsv(as.data.frame("No GO pathway found"), 
  gsea_go_file, 
  col_names = False)

}

if(nrow(gsea_kegg_result@result) > 0) {

  write_tsv(gsea_kegg_result@result, gsea_kegg_file)
  
  ## KEGG Dotplot
  pdf(file.path(out_dir, "gsea_kegg_dge_all_dotplot.pdf"), width = 8, height = 6)
  dotplot(gsea_kegg_result, showCategory = 20)
  dev.off()
  
  ## KEGG Ridgeplot
  pdf(file.path(out_dir, "gsea_kegg_dge_all_ridgeplot.pdf"), width = 8, height = 6)
  ridgeplot(gsea_kegg_result, showCategory = 20)
  dev.off()

  ## KEGG Top Pathway
  top_kegg_id <- gsea_kegg_result@result$ID[1]
  pdf(file.path(out_dir, "gsea_kegg_dge_all_top_pathway.pdf"), width = 8, height = 6)
  gseaplot2(gsea_kegg_result, geneSetID = top_kegg_id)
  dev.off()

} else {

  write_tsv(as.data.frame("No GO pathway found"), 
  gsea_kegg_file, 
  col_names = False)

}

##########################
# For up-regulated genes #
##########################
gene_list_filtered <- gene_list[gene_list > 0]
gene_list_filtered <- sort(gene_list_filtered, decreasing = TRUE)

gsea_go_result <- gsea_go(gene_list_filtered)
gsea_kegg_result <- gsea_kegg(gene_list_filtered)

gsea_go_file <- file.path(out_dir, "gsea_go_dge_up_results.tsv")
gsea_kegg_file <- file.path(out_dir, "gsea_kegg_dge_up_results.tsv")

if(nrow(gsea_go_result@result) > 0) {

  write_tsv(gsea_go_result@result, gsea_go_file)

} else {

  write_tsv(as.data.frame("No up GO pathway found"), 
  gsea_go_file, 
  col_names = False)

}

if(nrow(gsea_kegg_result@result) > 0) {

  write_tsv(gsea_kegg_result@result, gsea_kegg_file)

} else {

  write_tsv(as.data.frame("No up KEGG pathway found"), 
  gsea_kegg_file, 
  col_names = False)

}

############################
# For down-regulated genes #
############################
gene_list_filtered <- gene_list[gene_list < 0]
gene_list_filtered <- sort(gene_list_filtered, decreasing = FALSE)

gsea_go_result <- gsea_go(gene_list_filtered)
gsea_kegg_result <- gsea_kegg(gene_list_filtered)

if(nrow(gsea_go_result@result) > 0) {

  gsea_go_file <- file.path(out_dir, "gsea_go_dge_down_results.tsv")
  write_tsv(gsea_go_result@result, gsea_go_file)

} else {

  write_tsv(as.data.frame("No down GO pathway found"), 
  gsea_go_file, 
  col_names = False)

}

if(nrow(gsea_kegg_result@result) > 0) {

  gsea_kegg_file <- file.path(out_dir, "gsea_kegg_dge_down_results.tsv")
  write_tsv(gsea_kegg_result@result, gsea_kegg_file)

} else {

  write_tsv(as.data.frame("No down KEGG pathway found"), 
  gsea_kegg_file, 
  col_names = False)

}
