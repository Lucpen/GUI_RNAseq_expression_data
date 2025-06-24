# Module for filtering data based on annotation and selected accession IDs
# This module allows filtering of data based on the selected material and accession IDs.
# It provides two reactive data frames: one for all samples with the selected material
# and another for samples matching the selected accession IDs.

annotationDataFilteringServer <- function(id, df_filtered_annot, df_data, selected_acc_id, selected_control_type) {
    moduleServer(id, function(input, output, session) {

        # All samples with selected material (violin plot)
        df_filtered_data_ctrl <- reactive({
            req(df_filtered_annot(), df_data())
            annot <- df_filtered_annot()
            # Filter controls based on control_type
            if (selected_control_type() == "unaffected" && "affected" %in% names(annot)) {
                annot <- annot[!tolower(annot$affected) %in% "affected", , drop = FALSE]
            }
            if (is.null(selected_acc_id()) || length(selected_acc_id()) == 0) {
                df_filtered_annot_without_selected <- annot
            } else {
                df_filtered_annot_without_selected <- annot[!annot$InternalID %in% selected_acc_id(), , drop = FALSE]
            }
            sample_cols <- intersect(colnames(df_data()), as.character(df_filtered_annot_without_selected$InternalID))
            keep_cols <- c("ensembl_id", "hgnc_symbol", sample_cols)
            keep_cols <- intersect(keep_cols, colnames(df_data()))
            print("DEBUG: df_data colnames:")
            print(colnames(df_data()))
            print("DEBUG: InternalIDs in annotation:")
            print(as.character(df_filtered_annot_without_selected$InternalID))
            print("DEBUG: keep_cols:")
            print(keep_cols)
            if (length(sample_cols) == 0) {
                print("DEBUG: No sample columns matched, returning NULL")
                return(NULL)
            }
            result <- df_data()[, keep_cols, drop = FALSE]
            print("Filtered control data:")
            print(dim(result))
            print(head(result))
            result
        })
        
        # Only samples matching selected acc_id (dots)
        df_filtered_data_selected <- reactive({
            req(df_data())
            acc <- selected_acc_id()
            if (is.null(acc) || length(acc) == 0) {
                return(NULL)
            }
            sample_cols <- intersect(colnames(df_data()), as.character(acc))
            keep_cols <- c("ensembl_id", "hgnc_symbol", sample_cols)
            keep_cols <- intersect(keep_cols, colnames(df_data()))
            if (length(sample_cols) == 0) return(NULL)
            df_data()[, keep_cols, drop = FALSE]
        })
        
        return(list(
            df_filtered_data_ctrl = df_filtered_data_ctrl,
            df_filtered_data_selected = df_filtered_data_selected
        ))
    })
}

# DEBUG
# path_annot <- "/Users/luciapenaperez/workspace/CMMS/RNAseq/annotation_files/250610_samples_TPM_GUI.tsv"
#df_filtered_annot<-read.csv("/Users/luciapenaperez/workspace/CMMS/RNAseq/annotation_files/250610_samples_TPM_GUI.tsv")
# df_filtered_annot <- read.delim(path_annot)
#df_filtered_annot <- as.data.frame(df_filtered_annot)
#colnames(df_filtered_annot) <- trimws(colnames(df_filtered_annot))
#df_data <- readRDS("/Users/luciapenaperez/workspace/CMMS/RNAseq/data_files/250610_salmon_RNAseq_only_TPM.rds")
#sample_cols <- intersect(colnames(df_data), as.character(df_filtered_annot$IndividualID))