annotFilterUI <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(2,uiOutput(ns("material_ui"))),
    column(2,uiOutput(ns("control_type_ui"))),
    column(2, uiOutput(ns("family_ui"))),
    column(2, uiOutput(ns("individual_id_ui")))
    # REMOVE or comment out the next line:
    # column(3, uiOutput(ns("acc_id_ui")))
  )
}

annotFilterServer <- function(id, annot) {
  moduleServer(id, function(input, output, session) {

    # Reactive for filtered annotation by material
    filtered_annot <- reactive({
      annot_df <- annot()
      if (is.null(annot_df)) return(data.frame())
      if (!is.null(input$material) && input$material != "") {
        idx <- which(tolower(names(annot_df)) == "material")
        if (length(idx) > 0) {
          annot_df <- annot_df[annot_df[[idx[1]]] == input$material, , drop = FALSE]
        }
      }
    })

    output$material_ui <- renderUI({
      annot_df <- annot()
      selectizeInput(
        session$ns("material"),
        "Material",
        choices = if (is.null(annot_df)) "Loading..." else get_choices("Material", annot_df),
        multiple = FALSE,
        selected = input$material,
        options = list(placeholder = "Select Material")
      )
    })

    output$family_ui <- renderUI({
      annot_df <- annot()
      # Optionally filter by material if you want, but NOT by family
      if (!is.null(input$material) && input$material != "") {
        idx <- which(tolower(names(annot_df)) == "material")
        if (length(idx) > 0) {
          annot_df <- annot_df[annot_df[[idx[1]]] == input$material, , drop = FALSE]
        }
      }
      selectizeInput(
        session$ns("family"),
        "Family",
        choices = if (is.null(annot_df)) "Loading..." else get_choices("Family", annot_df),
        multiple = TRUE,
        selected = input$family,
        options = list(placeholder = "Select Family")
      )
    })

    # Update Individual choices based on selected family
    observeEvent(filtered_annot(), {
      valid_choices <- get_choices("IndividualID", filtered_annot())
      still_selected <- intersect(input$individual_id, valid_choices)
      updateSelectizeInput(
        session,
        "individual_id",
        choices = valid_choices,
        selected = still_selected
      )
    })

    output$individual_id_ui <- renderUI({
      annot_df <- filtered_annot()
      # Filter by selected families if any
      if (!is.null(input$family) && length(input$family) > 0) {
        idx <- which(tolower(names(annot_df)) == "family")
        if (length(idx) > 0) {
          annot_df <- annot_df[annot_df[[idx[1]]] %in% input$family, , drop = FALSE]
        }
      }
      selectizeInput(
        session$ns("individual_id"),
        "Individual",
        choices = if (nrow(annot_df) == 0) character(0) else sort(unique(annot_df$IndividualID)),
        multiple = TRUE,
        selected = input$individual_id,
        options = list(placeholder = "Select Individual")
      )
    })

    # Reactive for ACC numbers based on selected individuals
    acc_choices <- reactive({
      annot_df <- filtered_annot()
      if (!is.null(input$individual_id) && length(input$individual_id) > 0) {
        unique(annot_df$InternalID[annot_df$IndividualID %in% input$individual_id])
      } else {
        NULL
      }
    })

    output$control_type_ui <- renderUI({
      selectInput(
        session$ns("control_type"),
        "Use as controls:",
        choices = c("Unaffected" = "unaffected", "All samples" = "all"),
        selected = input$control_type %||% "unaffected"
      )
    })

    return(list(
      selected_material = reactive(input$material),
      selected_family = reactive(input$family),
      selected_individual_id = reactive(input$individual_id),
      selected_acc_id = acc_choices,
      selected_control_type = reactive(input$control_type),
      filtered_annot = reactive(filtered_annot())
    ))
  })
}