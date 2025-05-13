log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(GenomicFeatures)
library(txdbmaker)

# Replace with your GTF path
gtf_file <- snakemake@input[["gtf"]]
raw_counts <- read.table(snakemake@input[["counts"]], row.names = 1, header = TRUE)

txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")

exons_by_gene <- exonsBy(txdb, by = "gene")
gene_lengths <- sum(width(reduce(exons_by_gene)))
gene_lengths_kb <- gene_lengths / 1000

# Subset gene lengths to match raw_counts
common_genes <- intersect(rownames(raw_counts), names(gene_lengths))
raw_counts <- raw_counts[common_genes, ]
gene_lengths <- gene_lengths[common_genes]

rpk <- sweep(raw_counts, 1, gene_lengths / 1000, FUN = "/")

# Scaling factors (per sample)
scaling_factors <- colSums(rpk)

tpm <- sweep(rpk, 2, scaling_factors, FUN = "/") * 1e6
tpm$gene <- rownames(tpm)
tpm <- tpm[, c("gene", colnames(tpm)[1:(ncol(tpm)-1)])] # rearrange columns

write.table(tpm, file = snakemake@output[[1]], sep = "\t", quote = FALSE, row.names=FALSE)
