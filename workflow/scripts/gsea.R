# --- Load Required Packages ---
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db) # check for specific species when needed
  library(enrichplot)
  library(DOSE)
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(ggridges)
})

options(error = function() {
  traceback(2)
  quit(status = 1)
})

log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type="message")

out_dir <- snakemake@output[["output_dir"]]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# --- Read DESeq2 Results ---
df <- read_tsv(snakemake@input[["diffexpr"]])

# Expecting columns: gene, log2FoldChange
df <- df %>% dplyr::filter(!is.na(log2FoldChange))

# --- Convert Gene Symbols (ORFs) to Entrez IDs ---
gene_df <- bitr(df$gene, fromType = "ENSEMBL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)
print(head(gene_df))
# --- Merge Data and Rank Genes ---
merged_df <- inner_join(df, gene_df, by = c("gene" = "ENSEMBL"))
gene_list <- merged_df$log2FoldChange
names(gene_list) <- merged_df$ENTREZID
gene_list <- gene_list[!is.na(names(gene_list))]
dup_ids <- names(gene_list)[duplicated(names(gene_list))]
print(dup_ids)
gene_list <- gene_list[!duplicated(names(gene_list))]
gene_list <- sort(gene_list, decreasing = TRUE)

# --- Functions for pathway analysis ---
# --- (GO) ---
gsea_go <- function(gene_list, up_or_down) {
    gseGO(
        geneList = gene_list,
        OrgDb = org.Hs.eg.db,
        ont = "BP",
        minGSSize = 10,
        pvalueCutoff = 0.05,
        scoreType = up_or_down,
        verbose = FALSE
        )
}

# --- (KEGG) ---
gsea_kegg <- function(gene_list, up_or_down) {
    gseKEGG(
        geneList = gene_list,
        organism = 'hsa',  # KEGG code for human
        minGSSize = 10,
        pvalueCutoff = 0.05,
        scoreType = up_or_down,
        verbose = FALSE
        )

}


##############################
# For all disregulated genes #
##############################

gsea_go_result <- gsea_go(gene_list, "std")
gsea_kegg_result <- gsea_kegg(gene_list, "std")

gsea_go_file <- file.path(out_dir, "gsea_go_dge_all_results.tsv")
gsea_kegg_file <- file.path(out_dir, "gsea_kegg_dge_all_results.tsv")

if(nrow(gsea_go_result@result) > 0) {

  write_tsv(gsea_go_result@result, gsea_go_file)

  ## GO Dotplot
  pdf(file.path(out_dir, "gsea_go_dge_all_dotplot.pdf"), width = 8, height = 6)
  tryCatch({
    print(dotplot(gsea_go_result, showCategory = 20)) # explicitly force plot in the if loop
    }, finally = {
      dev.off()
      })
  
  ## GO Ridgeplot
  pdf(file.path(out_dir, "gsea_go_dge_all_ridgeplot.pdf"), width = 8, height = 6)
  tryCatch({
    print(ridgeplot(gsea_go_result, showCategory = 20))
    }, finally = {
    dev.off()
    })

  ## GO Top Pathway
  top_go_id <- gsea_go_result@result$ID[1]
  pdf(file.path(out_dir, "gsea_go_dge_all_top_pathway.pdf"), width = 8, height = 6)
  tryCatch({
    print(gseaplot(gsea_go_result, geneSetID = top_go_id))
    }, finally = {
    dev.off()
    })

} else {

  write_tsv(as.data.frame("No GO pathway found"), 
  gsea_go_file, 
  col_names = FALSE)

}

if(nrow(gsea_kegg_result@result) > 0) {

  write_tsv(gsea_kegg_result@result, gsea_kegg_file)
  
  ## KEGG Dotplot
  pdf(file.path(out_dir, "gsea_kegg_dge_all_dotplot.pdf"), width = 8, height = 6)
  tryCatch({
    print(dotplot(gsea_kegg_result, showCategory = 20))
    }, finally = {
    dev.off()
    })
  
  ## KEGG Ridgeplot
  pdf(file.path(out_dir, "gsea_kegg_dge_all_ridgeplot.pdf"), width = 8, height = 6)
  tryCatch({
    print(ridgeplot(gsea_kegg_result, showCategory = 20))
    }, finally = {
    dev.off()
    })

  ## KEGG Top Pathway
  top_kegg_id <- gsea_kegg_result@result$ID[1]
  pdf(file.path(out_dir, "gsea_kegg_dge_all_top_pathway.pdf"), width = 8, height = 6)
  tryCatch({
    print(gseaplot(gsea_kegg_result, geneSetID = top_kegg_id))
    }, finally = {
    dev.off()
    })

} else {

  write_tsv(as.data.frame("No GO pathway found"), 
  gsea_kegg_file, 
  col_names = FALSE)

}

##########################
# For up-regulated genes #
##########################

gsea_go_result <- gsea_go(gene_list, "pos")
gsea_kegg_result <- gsea_kegg(gene_list, "pos")

gsea_go_file <- file.path(out_dir, "gsea_go_dge_up_results.tsv")
gsea_kegg_file <- file.path(out_dir, "gsea_kegg_dge_up_results.tsv")

if(nrow(gsea_go_result@result) > 0) {

  write_tsv(gsea_go_result@result, gsea_go_file)

} else {

  write_tsv(as.data.frame("No up GO pathway found"), 
  gsea_go_file, 
  col_names = FALSE)

}

if(nrow(gsea_kegg_result@result) > 0) {

  write_tsv(gsea_kegg_result@result, gsea_kegg_file)

} else {

  write_tsv(as.data.frame("No up KEGG pathway found"), 
  gsea_kegg_file, 
  col_names = FALSE)

}

############################
# For down-regulated genes #
############################

gsea_go_result <- gsea_go(gene_list, "neg")
gsea_kegg_result <- gsea_kegg(gene_list, "neg")

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
  col_names = FALSE)

}
