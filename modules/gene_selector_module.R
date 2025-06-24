geneSearchUI <- function(id) {
  ns <- NS(id)
  fluidRow(
    # Left column: Ensembl ID input
    column(6,
      selectizeInput(
      ns("ensembl"),
      "Ensembl ID:",
      choices = NULL,
      options = list(
      placeholder = "Type Ensembl ID...",
      maxOptions  = 30
      )
      )
    ),
    # Right column: HGNC symbol input
    column(6,
      selectizeInput(
      ns("hgnc"),
      "HGNC Symbol:",
      choices = NULL,
      options = list(
      placeholder = "Type HGNC symbol...",
      maxOptions  = 30
      )
      )
    )
  )
}


geneSearchServer <- function(id, data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # When data() changes, populate dropdowns
    observeEvent(data(), {
      df <- data()
      updateSelectizeInput(session, "ensembl",
        choices = sort(unique(na.omit(df$ensembl_id))),
        server = TRUE)
      updateSelectizeInput(session, "hgnc",
        choices = sort(unique(na.omit(df$hgnc_symbol))),
        server = TRUE)
    })
    
    # Track last changed input
    last_changed <- reactiveVal("ensembl")
    
    observeEvent(input$ensembl, {
      req(data(), input$ensembl)
      last_changed("ensembl")
      
      df <- data()
      matched_hgnc <- unique(na.omit(df$hgnc_symbol[df$ensembl_id == input$ensembl]))
      
      if (length(matched_hgnc) > 0 && !is.na(matched_hgnc[1])) {
        if (input$hgnc != matched_hgnc[1]) {
          updateSelectizeInput(session, "hgnc", selected = matched_hgnc[1])
        }
      } else {
        updateSelectizeInput(session, "hgnc", selected = "")
      }
    }, ignoreInit = TRUE)
    
    observeEvent(input$hgnc, {
      req(data(), input$hgnc)
      last_changed("hgnc")
      
      df <- data()
      matched_ens <- unique(na.omit(df$ensembl_id[df$hgnc_symbol == input$hgnc]))
      
      if (length(matched_ens) > 0 && !is.na(matched_ens[1])) {
        if (input$ensembl != matched_ens[1]) {
          updateSelectizeInput(session, "ensembl", selected = matched_ens[1])
        }
      } else {
        updateSelectizeInput(session, "ensembl", selected = "")
      }
    }, ignoreInit = TRUE)
    
    # Output selected pair
    reactive({
      df <- data()
      if (last_changed() == "ensembl" &&
          !is.null(input$ensembl) && input$ensembl != "") {
        matched_hgnc <- unique(na.omit(df$hgnc_symbol[df$ensembl_id == input$ensembl]))
        return(list(ensembl = input$ensembl, hgnc = matched_hgnc[1]))
      }
      if (last_changed() == "hgnc" &&
          !is.null(input$hgnc) && input$hgnc != "") {
        matched_ens <- unique(na.omit(df$ensembl_id[df$hgnc_symbol == input$hgnc]))
        return(list(ensembl = matched_ens[1], hgnc = input$hgnc))
      }
      NULL
    })
  })
}
