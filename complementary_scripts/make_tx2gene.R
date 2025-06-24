# Save this as make_tx2gene.R and run it in your terminal or RStudio

library(rtracklayer)

# Path to your GTF file
gtf_file <- "/Users/luciapenaperez/workspace/CMMS/RNAseq/scripts/fb/gencode.v46.primary_assembly.annotation.gtf"
# Import GTF
gtf <- import(gtf_file)

# Extract transcript-to-gene mapping
tx2gene <- as.data.frame(mcols(gtf)[, c("transcript_id", "gene_id")])
tx2gene <- tx2gene[!is.na(tx2gene$transcript_id) & !is.na(tx2gene$gene_id), ]
colnames(tx2gene) <- c("TXNAME", "GENEID")
tx2gene$TXNAME <- sub("\\.[0-9]+$", "", tx2gene$TXNAME)

# Save as CSV
# write.csv(tx2gene, "tx2gene.csv", row.names = FALSE)

# Or save as RDS (recommended for R apps)
saveRDS(tx2gene, file = "tx2gene.rds")
