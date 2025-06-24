library(tximport)
library(rtracklayer)
library(biomaRt)

# Start by running this part in hasta
# Run from work dir containing all the quant.sf file
# You can also hav there the file with the paths to the files to export
main_path="data/250623_samples_TPM_blood_fibroblasts_and_muscle.txt"
gtf_file<-"references/v2.0/GRCh38/gencode.v46.primary_assembly.annotation.gtf"
date_="250623"

# Import GTF and create tx2gene mapping
gtf <- import(gtf_file)
tx2gene <- as.data.frame(mcols(gtf)[, c("transcript_id", "gene_id")])
tx2gene <- na.omit(tx2gene)
colnames(tx2gene) <- c("TXNAME", "GENEID")
tx2gene$TXNAME <- sub("\\.[0-9]+$", "", tx2gene$TXNAME)

# Prepare sample files
df_new_files <- read.csv(main_path, sep = "\t")
new_files <- df_new_files$salmon

# Import quantifications
txi <- tximport(new_files, type = "salmon", txIn = TRUE, tx2gene = tx2gene,
                countsFromAbundance = "no", ignoreTxVersion = TRUE)
for (assay in c("counts", "abundance", "length")) {
  colnames(txi[[assay]]) <- df_new_files$InternalID
}

# Prepare expression dataframe
df <- as.data.frame(txi$abundance)
df$ensembl_id <- sub("\\..*", "", rownames(df))

# Map Ensembl IDs to HGNC symbols
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
mapping <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = df$ensembl_id,
  mart = ensembl
)

df_expr <- merge(df, mapping, by.x = "ensembl_id", by.y = "ensembl_gene_id", all.x = TRUE)

# Save result
saveRDS(df, file = paste0(date_, "_salmon_RNAseq_only_TPM_muscle_fibroblast_blood.rds"))
