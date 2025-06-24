# Helper functions for the RNAseq GUI

# This function retrieves unique, sorted choices from a specified column in a filtered DataFrame.
# It returns "No column found" if the specified column does not exist in the DataFrame
get_choices <- function(col, filtered_df) {
      idx <- which(tolower(names(filtered_df)) == tolower(col))
      if (length(idx) == 0) return("No column found")
      sort(unique(na.omit(filtered_df[[idx[1]]])))
    }


# Helper function for all violin plots
make_violin_plot <- function(gene_row, annot, sel_long = NULL) {
  sample_cols <- setdiff(names(gene_row), c("ensembl_id", "hgnc_symbol"))
  long_df <- tidyr::pivot_longer(
    gene_row,
    cols = tidyselect::all_of(sample_cols),
    names_to = "sample",
    values_to = "expression"
  )
  print(paste("DEBUG: long_df rows:", nrow(long_df)))

  material <- NA
  if (!is.null(annot) && "InternalID" %in% names(annot) && "Material" %in% names(annot)) {
    material_row <- annot[annot$InternalID == sample_cols[1], "Material"]
    if (length(material_row) > 0) material <- material_row[1]
  }
  n_samples <- nrow(long_df)
  plot_title <- paste0(
    gene_row$hgnc_symbol, "\n(", gene_row$ensembl_id, ")\n",
    "Material: ", material, "\n",
    "Samples: ", n_samples
  )
  p <- ggplot(long_df, aes(x = "", y = expression)) +
    geom_violin(fill = "#A9CCE3", color = "black", trim = FALSE) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
    labs(
      title = plot_title,
      x = "Samples",
      y = "Expression (TPM)"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.text.x = element_blank(),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.text = element_text(size = 10)
    )
  if (!is.null(sel_long) && nrow(sel_long) > 0) {
    p <- p +
      geom_point(
        data = sel_long,
        aes(x = "", y = expression, color = IndividualID),
        size = 3,
        show.legend = TRUE
      )
  }
  p
}

# Helper function to get selected long format data
# This function retrieves the selected gene's data in long format
# and merges it with the annotation data if available.
# It returns NULL if the selected_df is NULL or if no matching ensembl_id is found
# in the selected_df.
get_selected_long <- function(selected_df, annot, ensembl_id) {
  if (is.null(selected_df)) return(NULL)
  row <- selected_df[selected_df$ensembl_id == ensembl_id, , drop = FALSE]
  if (nrow(row) == 0) return(NULL)
  sample_cols <- setdiff(names(row), c("ensembl_id", "hgnc_symbol"))
  if (length(sample_cols) == 0) return(NULL)
  sel_long <- tidyr::pivot_longer(
    row,
    cols = tidyselect::all_of(sample_cols),
    names_to = "sample",
    values_to = "expression"
  )
  if (!is.null(annot) && all(c("InternalID", "IndividualID", "Material") %in% names(annot))) {
    sel_long <- merge(
      sel_long,
      annot[, c("InternalID", "IndividualID", "Material")],
      by.x = "sample",
      by.y = "InternalID",
      all.x = TRUE
    )
  }
  sel_long
}

check_data_to_annotation_consistency <- function(df, annot) {
  colnames_df <- colnames(df)
  colnames_df_internalID <- grep("ACC", colnames_df, value = TRUE, ignore.case = TRUE)
  internalID_annot <- annot$InternalID
  overlap <- intersect(colnames_df_internalID, internalID_annot)
  internalid_not_in_data <- NULL
  individualid_not_in_data <- NULL
  if (length(overlap) == 0) {
    message_for_user <- "No matching InternalID found between data and annotation files."
  } else if (length(overlap) < length(colnames_df_internalID)) {
    not_mathing <- setdiff(colnames_df_internalID, internalID_annot)
    message_for_user <- paste0( "Some InternalIDs in the data do not match those in the annotation file: ",
            paste(not_mathing, collapse = ", "))
  } else if (length(overlap) < length(internalID_annot)) {
    not_mathing <- setdiff(internalID_annot, colnames_df_internalID)
    internalid_not_in_data <- annot$InternalID[annot$InternalID %in% not_mathing]
    individualid_not_in_data <- annot$IndividualID[annot$InternalID %in% not_mathing]
    message_for_user <- paste0( 
      "Some InternalIDs in the annotation file do not match those in the data: ",
        paste(not_mathing, collapse = ", "))
  } else {
    message_for_user <- "All InternalIDs in the data match those in the annotation file."
  }
  return(list(
    message = message_for_user,
    overlap = overlap,
    individualid_not_in_data = individualid_not_in_data,
    internalid_not_in_data = internalid_not_in_data
  ))
}

# Helper function to filter by gene
filter_by_gene <- function(df, sel) {
  if (is.null(df) || is.null(sel)) return(NULL)
  if (!is.null(sel$ensembl) && sel$ensembl != "") {
    df[df$ensembl_id == sel$ensembl, , drop = FALSE]
  } else if (!is.null(sel$hgnc) && sel$hgnc != "") {
    df[df$hgnc_symbol == sel$hgnc, , drop = FALSE]
  } else {
    NULL
  }
}

add_samples_to_data <- function(annotation_data, loaded_data, internalid_not_in_data, tx2gene_file, import_level = "gene") {
  library(tximport)

  # Load precomputed tx2gene (supports .rds or .csv)
  if (grepl("\\.rds$", tx2gene_file, ignore.case = TRUE)) {
    tx2gene <- readRDS(tx2gene_file)
  } else if (grepl("\\.csv$", tx2gene_file, ignore.case = TRUE)) {
    tx2gene <- read.csv(tx2gene_file, stringsAsFactors = FALSE)
  } else {
    stop("tx2gene_file must be a .rds or .csv file")
  }

  # List with all new samples
  df_new_files <- annotation_data %>%
    dplyr::filter(InternalID %in% internalid_not_in_data)
  new_files <- df_new_files$salmon

  # Set tximport arguments and run tximport based on import_level
  if (import_level == "gene") {
    tximport_args <- list(type = "salmon", txIn = TRUE, tx2gene = tx2gene, ignoreTxVersion = TRUE)
    tximport_args_TPM <- c(tximport_args, countsFromAbundance = "scaledTPM")
    txi_new <- do.call(tximport, c(list(files = new_files), tximport_args_TPM))
    # Remove the version suffix from Ensembl IDs
    # This is necessary to match the loaded_data Ensembl IDs and to select them correctly afterwards
    rownames(txi_new$abundance) <- gsub("\\..*", "", rownames(txi_new$abundance))
  } else if (import_level == "transcript") {
    tximport_args <- list(type = "salmon", txIn = TRUE, txOut = TRUE, ignoreTxVersion = TRUE)
    tximport_args_TPM <- c(tximport_args, countsFromAbundance = "scaledTPM")
    txi_new <- do.call(tximport, c(list(files = new_files), tximport_args_TPM))

  } else {
    stop("Invalid import_level: must be 'gene' or 'transcript'")
  }

  colnames(txi_new$abundance) <- df_new_files$InternalID
  txi_new$abundance <- cbind(ensembl_id = rownames(txi_new$abundance), as.data.frame(txi_new$abundance))
  # Merge new samples with existing data (by rownames)
  complete_data <- merge(
    txi_new$abundance,
    loaded_data,
    by = "ensembl_id",
    all = TRUE
  ) %>% 
  dplyr::select(
      ensembl_id, hgnc_symbol, dplyr::everything()
    )

  # Check if the number of rows matches to make sure same reference was used
  if (nrow(complete_data) != nrow(loaded_data)) {
    stop("The length of the new data does not match the length of the loaded data.")
  }

  return(complete_data)
}
