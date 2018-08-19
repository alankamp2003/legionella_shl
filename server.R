library(shiny)
library(shinyjs)
library(readxl)
library(xlsx)
library(dplyr)
library(DT)
source("legionellaFuncs.R")

init_data <- NULL

# enables the download button with the passed key if the passed 
# data (a list) is not null and not empty; otherwise disables the button
toggleDLButton <- function(key, data) {
  toggleState(key, !is.null(data) && length(data) > 0)
}

# validates the data in the input fields
validateInput <- function(input) {
  if (is.null(input$file1)) {
    "Please choose a legionella data file"
  } else if (!is.numeric(input$coverage)) {
    "Please enter a valid number as Coverage"
  } else if (!is.numeric(input$identity)) {
    "Please enter a valid number as Identity"
  } else {
     NULL
  }
}

# converts the passed vectors elements into ones to be shown in a select widget;
# does so by gets unique elements from the passed vector, sorting them and addin
# the "All" option to them
getSelectOpts <- function(options) {
  options <- sort(unique(options))
  options <- c("All" = "", options)
}

shinyServer(function(input, output, session) {
  observe({
    inFile <- input$file1
    if (is.null(inFile)) {
      toggleDLButton('download_all', NULL)
      toggleDLButton('download_genes', NULL)
      return(NULL)
    }
    # initialize the data structure containing all data only once
    if (is.null(init_data)) {
      # Create a Progress object
      progress <- shiny::Progress$new()
      progress$set(message = "Loading data. Please wait...", value = 0)
      # Close the progress when this reactive exits (even if there's an error)
      on.exit(progress$close())
      
      # Create a callback function to update progress.
      # Each time this is called:
      # - If `value` is NULL, it will move the progress bar 1/5 of the remaining
      #   distance. If non-NULL, it will set the progress to that value.
      # - It also accepts optional detail text.
      updateProgress <- function(value = NULL, detail = NULL) {
        if (is.null(value)) {
          value <- progress$getValue()
          value <- value + (progress$getMax() - value) / 5
        }
        progress$set(value = value, detail = detail)
      }
      init_data <<- getInitialData(inFile$datapath, updateProgress)
      
      # load the select(drop down) widgets from the data
      species <- getSelectOpts(init_data$Species)
      updateSelectInput(session, "species",
                        choices = species)
      serogroups <- getSelectOpts(init_data$Serogroup)
      updateSelectInput(session, "serogroup",
                        choices = serogroups)
      seq_types <- getSelectOpts(init_data$`Sequence Type`)
      updateSelectInput(session, "seq_type",
                        choices = seq_types)
      genes <- getSelectOpts(init_data$GENE)
      args <- read_excel('args.xlsx')
      sel_genes <- args[[1, 'Genes']]
      if (!is.na(sel_genes))
        sel_genes <- strsplit(sel_genes, ',')[[1]]
      else
        sel_genes <- c()
      updateSelectInput(session, "gene",
                        choices = genes, selected = sel_genes)
    }
  })
  
  tables <- eventReactive(input$gen_tables, {
    # Create a Progress object
    progress <- shiny::Progress$new()
    progress$set(message = "Generating output...", value = 0)
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())
    
    # Create a callback function to update progress.
    # Each time this is called:
    # - If `value` is NULL, it will move the progress bar 1/5 of the remaining
    #   distance. If non-NULL, it will set the progress to that value.
    # - It also accepts optional detail text.
    updateProgress <- function(value = NULL, detail = NULL) {
      if (is.null(value)) {
        value <- progress$getValue()
        value <- value + (progress$getMax() - value) / 5
      }
      progress$set(value = value, detail = detail)
    }
    
    # pass a copy of the original data so that filtering doesn't change it
    curr_data <- init_data[,]
    getFilteredData(curr_data, species = input$species, serogroup = input$serogroup, seq_type = input$seq_type,
                    coverage = input$coverage, identity = input$identity, 
                    genes = input$gene, updateProgress = updateProgress)
  })
  
  output$all_fields <- renderDataTable({
    validate(validateInput(input))
    fields <- tables()$all_fields
    toggleDLButton('download_all', fields)
    datatable(fields, options = list(pageLength = 10), rownames= FALSE)
  })

  output$gene_fields <- renderDataTable({
    validate(validateInput(input))
    fields <- tables()$gene_fields
    toggleDLButton('download_genes', fields)
    datatable(fields, options = list(pageLength = 10), rownames= FALSE)
  })
  
  # Downloadable xslx of all fields
  output$download_all <- downloadHandler(
    filename = function() {
      paste("all_fields", ".xlsx", sep = "")
    },
    content = function(file) {
      df <- data.frame(tables()$all_fields)
      write.xlsx(df, file, row.names = FALSE)
    }
  )
  
  # # Downloadable xslx of gene fields
  output$download_genes <- downloadHandler(
    filename = function() {
      paste("gene_fields", ".xlsx", sep = "")
    },
    content = function(file) {
      df <- data.frame(tables()$gene_fields)
      write.xlsx(df, file, row.names = FALSE)
    },
    contentType = "application/vnd.ms-excel"
  )
})
