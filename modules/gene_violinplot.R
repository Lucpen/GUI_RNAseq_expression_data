
# Minimal violin plot module (for housekeeper tab, plot only)
violin_plot_UI <- function(id) {
  ns <- NS(id)
  plotOutput(ns("violin_plot"))
}

violin_plot_Server <- function(id, gene_row, annot, sel_long = NULL) {
  moduleServer(id, function(input, output, session) {
    output$violin_plot <- renderPlot({
      make_violin_plot(gene_row, annot, sel_long)
    })
  })
}

# Full-featured gene plot module (for main gene expression tab)
gene_plot_UI <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(
      width = 4,
      plotOutput(ns("violin_plot"))
    ),
    column(
      width = 3,
      tableOutput(ns("ctrl_stats_table")),
      br(),
      tableOutput(ns("selected_gene_long_table")),
      br(),
      violin_downloads_UI(ns("downloads"))
    )
  )
}

gene_plot_Server <- function(id, ctrl_data, selected_data, filtered_annot) {
  moduleServer(id, function(input, output, session) {
    sel_long_data <- reactive({
      sel_df <- selected_data()
      annot <- filtered_annot()
      # Use the ensembl_id from sel_df if available
      ensembl_id <- if (!is.null(sel_df) && "ensembl_id" %in% names(sel_df)) sel_df$ensembl_id[1] else NULL
      get_selected_long(sel_df, annot, ensembl_id)
    })

    violin_plot_obj <- reactive({
      df <- ctrl_data()
      req(df)
      gene_row <- df
      validate(
        need(nrow(gene_row) == 1, "Please select a single gene.")
      )
      annot <- filtered_annot()
      sel_long <- sel_long_data()
      make_violin_plot(gene_row, annot, sel_long)
    })

    output$violin_plot <- renderPlot({
      violin_plot_obj()
    })

    output$selected_gene_long_table <- renderTable({
      validate(
        need(!is.null(sel_long_data()), "Select a sample to see details.")
      )
      sel_long_data() %>% 
        dplyr::select(hgnc_symbol, Material, IndividualID, expression) %>%
        dplyr::arrange(IndividualID) %>% 
        dplyr::rename("Expression(TPM)" = "expression",
        "Hgnc symbol" = "hgnc_symbol")
    })

    ctrl_stats <- reactive({
      df <- ctrl_data()
      req(df)
      sample_cols <- setdiff(names(df), c("ensembl_id", "hgnc_symbol"))
      if (length(sample_cols) == 0) return(NULL)
      values <- as.numeric(as.matrix(df[, sample_cols, drop = FALSE]))
      data.frame(
        hgnc_symbol = df$hgnc_symbol[1],
        ensembl_id = df$ensembl_id[1],
        Mean = mean(values, na.rm = TRUE),
        Median = median(values, na.rm = TRUE),
        SD = sd(values, na.rm = TRUE)
      )
    })

    output$ctrl_stats_table <- renderTable({
      ctrl_stats() %>%
      dplyr::rename("Ensembl ID" = "ensembl_id",
        "Hgnc symbol" = "hgnc_symbol")
    }, rownames = FALSE)

    violin_downloads_Server(
      "downloads",
      plot_obj = violin_plot_obj,
      gene_row = reactive({ ctrl_data() }),
      annot = reactive({ filtered_annot() })
    )
  })
}

violin_downloads_UI <- function(id) {
  ns <- NS(id)
  div(
    style = "display: flex; gap: 10px; margin-top: 10px; width: 100%;",
    downloadButton(ns("download_plot_png"), "Download plot PNG", style = "width: 100%;"),
    downloadButton(ns("download_plot_svg"), "Download plot SVG", style = "width: 100%;")
  )
}

violin_downloads_Server <- function(id, plot_obj, gene_row, annot) {
  moduleServer(id, function(input, output, session) {
    output$download_plot_png <- downloadHandler(
      filename = function() {
        gene_row_val <- gene_row()
        annot_val <- annot()
        symbol <- if (!is.null(gene_row_val) && "hgnc_symbol" %in% names(gene_row_val)) gene_row_val$hgnc_symbol[1] else "gene"
        material <- if (!is.null(annot_val) && "Material" %in% names(annot_val)) annot_val$Material[1] else "material"
        paste0("violin_plot_", symbol, "_", material, "_", Sys.Date(), ".png")
      },
      content = function(file) {
        ggsave(file, plot = plot_obj(), width = 5, height = 6, dpi = 300, bg = "transparent")
      }
    )
    output$download_plot_svg <- downloadHandler(
      filename = function() {
        gene_row_val <- gene_row()
        annot_val <- annot()
        symbol <- if (!is.null(gene_row_val) && "hgnc_symbol" %in% names(gene_row_val)) gene_row_val$hgnc_symbol[1] else "gene"
        material <- if (!is.null(annot_val) && "Material" %in% names(annot_val)) annot_val$Material[1] else "material"
        paste0("violin_plot_", symbol, "_", material, "_", Sys.Date(), ".svg")
      },
      content = function(file) {
        svg(filename = file, width = 5, height = 6, bg = "transparent")
        print(plot_obj())
        dev.off()
      }
    )
  })
}
