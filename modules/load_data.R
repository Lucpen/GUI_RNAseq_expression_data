# Module to load data and annotation files

buttons_tab1_UI <- function(id) {
  ns <- NS(id)
  tagList(
    fileInput(inputId = ns("load_db"),
              label = "Load RDS file",
              accept = ".rds"),
    fileInput(inputId = ns("load_annot"),
              label = "Load annotation file",
              accept = c(".tsv", ".csv", ".txt", ".xlsx", ".xls")), # Accepts xcel, tsv and csv files
    uiOutput(ns("sample_info")),
    uiOutput(ns("buttons_upload")),
    uiOutput(ns("save_data_ui")) # Render download button conditionally
  )
}

buttons_tab1_Server <- function(id) {
  moduleServer(id, function(input, output, session) {
    # Reactive value to store loaded data
    loaded_data <- reactiveVal(NULL)
    annotation_data <- reactiveVal(NULL)
    data_updated <- reactiveVal(FALSE) # Track if data was updated

    # Load RData file
    observeEvent(input$load_db, {
      print("File upload observer fired!")
      req(input$load_db)
      df <- readRDS(input$load_db$datapath)
      df <- as.data.frame(df)  # Ensure it's a plain data frame
      colnames(df) <- trimws(colnames(df))  # Remove whitespace from column names
      loaded_data(df)
      data_updated(FALSE) # Reset flag on new load
    })

    # Load annotation file (.tsv, .csv, .txt, .xlsx, .xls)
    observeEvent(input$load_annot, {
      req(input$load_annot)
      ext <- tools::file_ext(input$load_annot$name)
      if (ext %in% c("xlsx", "xls")) {
        annot <- readxl::read_excel(input$load_annot$datapath)
        annot <- as.data.frame(annot)
      } else {
        annot <- read.delim(input$load_annot$datapath)
        annot <- as.data.frame(annot)
      }
      colnames(annot) <- trimws(colnames(annot))
      
      # Ensure 'affected' column exists and is standardized
      if (!"affected" %in% tolower(colnames(annot))) {
        annot$affected <- "affected"
      }
      # Standardize column name (in case it's 'Affected' or similar)
      affected_col <- which(tolower(colnames(annot)) == "affected")
      colnames(annot)[affected_col] <- "affected"
      annot$affected <- tolower(trimws(annot$affected))
      
      # Set to "unaffected" if InternalID contains "-II-"
      if ("IndividualID" %in% colnames(annot)) {
        annot$affected[grepl("-II-", annot$IndividualID)] <- "unaffected"
      }
      
      annotation_data(annot)
      data_updated(FALSE) # Reset flag on new load
    })

    data_check <- reactive({
      req(annotation_data(), loaded_data())
      check_data_to_annotation_consistency(
        loaded_data(),
        annotation_data()
      )
    })

    output$sample_info <- renderUI({
      req(annotation_data(), loaded_data(), data_check())
      HTML(paste( data_check()$message, "<br>",
        "A total of ", length(data_check()$overlap), "samples properly loaded"
      ))
    })

    output$buttons_upload <- renderUI({
      req(loaded_data(), annotation_data(), data_check())
      if (!is.null(data_check()$individualid_not_in_data)) {
        tagList(
          HTML(paste(
            "The following samples were not found in the data: ",
            paste(data_check()$individualid_not_in_data, collapse = ", "), "<br>",
            "Make sure the path given in the salmon files is correct and press the load button to add them to the data."
          )),
          actionButton(session$ns("load_data"), "Load missing samples")
        )
      }
    })

    observeEvent(input$load_data, {
      req(annotation_data(), loaded_data(), data_check())
      if (!is.null(data_check()$internalid_not_in_data)) {
        tx2gene_file <- "data/tx2gene.rds"

        withProgress(message = "Adding missing samples...", value = 0, {
          incProgress(0.2, detail = "Importing new samples...")
          updated_data <- add_samples_to_data(
            annotation_data = annotation_data(),
            loaded_data = loaded_data(),
            internalid_not_in_data = data_check()$internalid_not_in_data,
            tx2gene_file = tx2gene_file,
            import_level = "gene"
          )
          print(colnames(updated_data))  # Print column names for verification
          incProgress(0.8, detail = "Updating data...")
          loaded_data(updated_data)
        })

        showNotification("Missing samples added to data!", type = "message")
        data_updated(TRUE) # Set flag to TRUE after update
      }
    }, ignoreInit = TRUE)
    
    # Download handler for saving data as RDS
    output$save_data <- downloadHandler(
      filename = function() {
        paste0(Sys.Date(),"_salmon_RNAseq_only_TPM.rds")
      },
      content = function(file) {
        req(loaded_data())
        saveRDS(loaded_data(), file)
      }
    )

    # Conditionally render the download button
    output$save_data_ui <- renderUI({
      if (isTRUE(data_updated())) {
        downloadButton(session$ns("save_data"), "Save updated data as RDS")
      }
    })

    # Return the reactive values
    return(list(
      loaded_data = loaded_data,
      annotation_data = annotation_data
    ))
  })
}
