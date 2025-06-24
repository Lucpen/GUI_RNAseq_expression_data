housekeeperUI <- function(id) {
  ns <- NS(id)
  tagList(
    annotFilterUI(ns("housekeeper_filters")),
    fluidRow(
      column(4,
        plotOutput(ns("hk1_violin_plot"), height = "250px"),
        violin_downloads_UI(ns("hk1_dl"))
      ),
      column(4,
        plotOutput(ns("hk2_violin_plot"), height = "250px"),
        violin_downloads_UI(ns("hk2_dl"))
      )
    ),
    fluidRow(
      column(4,
        plotOutput(ns("hk3_violin_plot"), height = "250px"),
        violin_downloads_UI(ns("hk3_dl"))
      ),
      column(4,
        plotOutput(ns("hk4_violin_plot"), height = "250px"),
        violin_downloads_UI(ns("hk4_dl"))
      )
    )
  )
}

housekeeperServer <- function(id, filtered_data, gene_filters) {
  moduleServer(id, function(input, output, session) {
    hk_genes <- c("ENSG00000075624", "ENSG00000156508", "ENSG00000089157", "ENSG00000111640")
    ns <- session$ns

    lapply(seq_along(hk_genes), function(i) {
      local({
        idx <- i
        gene_id <- hk_genes[idx]
        sel <- list(ensembl = gene_id, hgnc = NULL)
        plot_id <- paste0("hk", idx, "_violin_plot")
        dl_id <- paste0("hk", idx, "_dl")

        gene_row <- reactive({
          filter_by_gene(filtered_data$df_filtered_data_ctrl(), sel)
        })
        sel_long <- reactive({
          get_selected_long(
            filtered_data$df_filtered_data_selected(),
            gene_filters$filtered_annot(),
            gene_id
          )
        })
        annot <- reactive({ gene_filters$filtered_annot() })

        output[[plot_id]] <- renderPlot({
          make_violin_plot(gene_row(), annot(), sel_long())
        })

        violin_downloads_Server(
          dl_id,
          plot_obj = reactive({ make_violin_plot(gene_row(), annot(), sel_long()) }),
          gene_row = gene_row,
          annot = annot
        )
      })
    })
  })
}