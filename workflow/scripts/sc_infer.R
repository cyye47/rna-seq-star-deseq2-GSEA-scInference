log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(granulator)

# reference dataset PBMC from 
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107011
# representing bulk PBMC expression levels from 13 individuals
# as well as 29 isolated immune cell types from 4 individuals

### this part run only once to generate an immune cell transcriptome reference
### they are based on mean expression of sorted cell populations from 4 individuals

# library(dplyr)
# library(tidyr)
# reference_data <- read.table("/Users/chaoyangye/Documents/Consulting/BridgeInfomatics/resources/hg38_chr22/GSE107011_Processed_data_TPM.txt",
#                              check.names = FALSE)
# 
# # process to get average of isolated immune cell types from the 4 individuals
# # this will be the true reference
# # first remove PBMC
# reference_data <- reference_data[, grep("PBMC", colnames(reference_data), invert = TRUE)]
# # then calculate average by 29 cell types
# # B_Ex, B_naive, B_NSM, B_SM, Basophils, C_mono
# # CD4_naive, CD4_TE, CD8_CM, CD8_EM, CD8_naive, CD8_TE
# # I_mono, MAIT, mDC, NC_mono, Neutrophils, NK, pDC, Plasmablasts
# # Progenitor, TFH, Th1, Th1/Th17, TH17, Th2, Treg, VD2-, VD2+
# 
# # Create a lookup table of column names and their group (based on suffix)
# column_groups <- tibble(
#   colname = colnames(reference_data),
#   group = sub(".*_(B_Ex|B_naive|B_NSM|B_SM|Basophils|C_mono|CD4_naive|CD4_TE|CD8_CM|CD8_EM|CD8_naive|CD8_TE|I_mono|MAIT|mDC|NC_mono|Neutrophils|NK|pDC|Plasmablasts|Progenitor|TFH|Th1|Th1/Th17|Th17|Th2|Treg|VD2\\-|VD2\\+|)$", 
#               "\\1", 
#               colnames(reference_data))
# )
# 
# # Compute row means by group
# group_rowmeans_df <- column_groups %>%
#   group_by(group) %>%
#   summarise(rowmean = list(rowMeans(reference_data[, colname]))) %>%
#   unnest_wider(rowmean)
#
# colnames(group_rowmeans_df)[1] <- "gene"
# to_write <- data.frame(t(group_rowmeans_df))
# colnames(to_write) <- to_write[1, ]
# to_write <- to_write[2:nrow(to_write), ]
# to_write$gene <- rownames(to_write)
# to_write <- to_write[, c("gene", colnames(to_write))[1:29]]
# write.csv(to_write, quote=FALSE, row.names=FALSE,
#           file = "/Users/chaoyangye/Documents/Consulting/BridgeInfomatics/resources/hg38_chr22/mean_expr_im_types.csv")

reference_mean <- read.csv(snakemake@input[["pbmc_reference"]],
                          header = TRUE, check.names = FALSE)
rownames(reference_mean) <- reference_mean$gene
reference_mean <- as.matrix(reference_mean[, 2:ncol(reference_mean)])

ref_signature <- list(ref = reference_mean)

# read bulkRNAseq data input in TPM count with ensembl ID
bulkRNAseq <- read.table(snakemake@input[["tpm"]],
                         header=TRUE)
rownames(bulkRNAseq) <- bulkRNAseq$gene
bulkRNAseq <- as.matrix(bulkRNAseq[, 2:ncol(bulkRNAseq)])

# deconvolute input data using all methods and reference profile matrix
methods <- c('ols','nnls','qprog','rls','svr', 'dtangle', 'qprogwc')
decon <- deconvolute(bulkRNAseq, ref_signature, methods)

# plot cell type proportions across models
# this helps to find concensus and model performance
pdf(snakemake@output[["model_performance"]], width = 8, height = 6)
tryCatch({
    plot_deconvolute(deconvoluted = decon, scale = TRUE, labels = FALSE)
    }, finally = {
    dev.off()
    })
# check for consistent performance across models

# correlation analysis comparing methods
correl <- correlate(deconvoluted = decon)
# rank the best performing method, which correlates the best with other methods
correl$rank
best_method <- correl$rank[1, "method"]
# option to choose an alternative method
# best_method <- "svr"

# plot cell type proportions for the best method
pdf(snakemake@output[["proportion_fig"]], width = 8, height = 6)
tryCatch({
    plot_proportions(deconvoluted = decon, 
                 method = best_method, 
                 signature = names(ref_signature))
    }, finally = {
    dev.off()
    })

# write cell proportions with the best method
best_model <- paste(best_method, names(ref_signature), sep="_")
decon$proportions[[best_model]]

write.table(decon$proportions[[best_model]], 
      file = snakemake@output[["proportion_table"]], 
      sep = "\t", 
      quote = FALSE)
