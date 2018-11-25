library(shiny)
library(shinyjs)
library(readxl)
library(xlsx)
library(dplyr)
library(DT)
source("legionellaFuncs.R")

init_data <- NULL
meta <- NULL

# enables the download button with the passed key if the passed 
# data (a list) is not null, not empty and none of the vectors 
# in it are empty; otherwise disables the button
toggleDLButton <- function(key, data) {
  togl <- !is.null(data) && length(data) > 0
  if (togl) {
    k <- sapply(data, function(v) length(v) > 0)
    togl <- sum(k) > 0
  }
  toggleState(key, togl)
}

# validates the data in the input fields
validateInput <- function(input) {
  if (is.null(input$file1)) {
    "Please choose a virulence data file"
  } else if (!is.numeric(input$coverage)) {
    "Please enter a valid number as Seq Total % Coverage"
  } else if (!is.numeric(input$identity)) {
    "Please enter a valid number as Seq Weighted % Identity"
  } else {
    NULL
  }
}

# converts the passed vectors elements into ones to be shown in a select widget;
# does so by getting unique elements from the passed vector, sorting them and adding
# the "All" option to them
getSelectOpts <- function(options) {
  options <- sort(unique(options))
  options <- c("All" = "", options)
}

# Stores the plots contained in the passed list as plotOutput objects so that they can be displayed
# later in "Visualization" tab; "prefix" is used so that plots can be found at the time of display
storePlots <- function(plot_list, prefix) {
  plot_output_list <- lapply(1:length(plot_list), function(i) {
    plotname <- paste(prefix, i, sep="")
    plotOutput(plotname, height = 280, width = 300, inline = TRUE)
  })
  do.call(tagList, plot_output_list)
}

# Displays the plots contained in the passed list in "Visualization" tab; "prefix" is used to 
# find specific plots that are stored in the passed output already 
displayPlots <- function(plot_list, prefix, output) {
  #print(paste("plots", length(plot_list)))
  if (length(plot_list) > 0) {
    for (i in 1:length(plot_list)) {
      # Need local so that each item gets its own number. Without it, the value
      # of i in the renderPlot() will be the same across all instances, because
      # of when the expression is evaluated.
      local({
        my_i <- i
        plotname <- paste(prefix, my_i, sep="")
        output[[plotname]] <- renderPlot(expr = {
          plot_list[[my_i]]
        }, height = 200, width = 300)
      })
    }
  }
}

shinyServer(function(input, output, session) {
  observe({
    inFile <- input$file1
    if (is.null(inFile)) {
      toggleDLButton('download_all', NULL)
      toggleDLButton('download_genes', NULL)
      toggleDLButton('download_vir_hits', NULL)
      toggleDLButton('download_stats_ana', NULL)
      #toggleDLButton('download_vis', NULL)
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
      return_list <- getInitialData(inFile$datapath, updateProgress)
      init_data <<- return_list$fields
      meta <<- return_list$meta
      
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
      # the commented out code reads the xlsx file; if the file contains a list of genes,
      # sets them as the genes selected by default in the dropdown for genes; the code 
      # can be uncommented later if need be
      # args <- read_excel('args.xlsx')
      # sel_genes <- args[[1, 'Genes']]
      # if (!is.na(sel_genes))
      #   sel_genes <- strsplit(sel_genes, ',')[[1]]
      # else
      #   sel_genes <- c()
      # updateSelectInput(session, "gene",
      #                   choices = genes, selected = sel_genes)
      updateSelectInput(session, "gene", choices = genes)
    }
  })
  
  # the following observeEvents manage the checkboxes for species, serogroup and sequence type, in that order;
  # the order also signifies parent child relationship, from left to right; if a child's checkbox is checked,
  # its parent's checkbox is automatically checked; if a parent's checkbox is unchecked, its child's checkbox is
  # is automatically unchecked
  observeEvent(input$vir_spec, {
    if (!input$vir_spec && input$vir_sero) {
      updateCheckboxInput(session, "vir_sero", value = FALSE)
    }
  })
  
  observeEvent(input$vir_sero, {
    if (input$vir_sero && !input$vir_spec) {
      updateCheckboxInput(session, "vir_spec", value = TRUE)
    } else if (!input$vir_sero && input$vir_seq_type) {
      updateCheckboxInput(session, "vir_seq_type", value = FALSE)
    }
  })
  
  observeEvent(input$vir_seq_type, {
    if (input$vir_seq_type && !input$vir_sero) {
      updateCheckboxInput(session, "vir_sero", value = TRUE)
    }
  })
  
  tables <- eventReactive(input$gen_output, {
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
    output <- getFilteredData(curr_data, meta, species = input$species, serogroup = input$serogroup, seq_type = input$seq_type,
                              coverage = input$coverage, identity = input$identity, 
                              genes = input$gene, show_spec = input$vir_spec, 
                              show_sero = input$vir_sero, show_seq_type = input$vir_seq_type, 
                              updateProgress = updateProgress)
    output
  })
  
  output$all_fields <- renderDataTable({
    toggleDLButton('download_all', NULL)
    validate(validateInput(input))
    fields <- tables()$all_fields
    toggleDLButton('download_all', fields)
    datatable(fields, options = list(pageLength = 10), rownames= FALSE)
  })
  
  output$gene_fields <- renderDataTable({
    toggleDLButton('download_genes', NULL)
    validate(validateInput(input))
    fields <- tables()$gene_fields
    toggleDLButton('download_genes', fields)
    datatable(fields, options = list(pageLength = 10), rownames= FALSE)
  })
  
  output$virulence_hits <- renderDataTable({
    toggleDLButton('download_vir_hits', NULL)
    validate(validateInput(input))
    fields <- tables()$virulence_hits
    toggleDLButton('download_vir_hits', fields)
    datatable(fields, options = list(pageLength = 10), rownames= FALSE)
  })
  
  output$stats_analysis <- renderDataTable({
    toggleDLButton('download_stats_ana', NULL)
    validate(validateInput(input))
    fields <- tables()$stats_analysis
    toggleDLButton('download_stats_ana', fields)
    datatable(fields, options = list(pageLength = 10), rownames= FALSE)
  })
  
  # output$plot1 <- renderPlot({
  #   tables()$plots$plot1
  # })
  
  output$spec_plots <- renderUI({
    storePlots(tables()$plots$spec_plots, "spec_plot")
  })
  
  output$sero_plots <- renderUI({
    storePlots(tables()$plots$sero_plots, "sero_plot")
  })
  
  output$seq_type_plots <- renderUI({
    storePlots(tables()$plots$seq_type_plots, "seq_type_plot")
  })
  
  # Downloadable xslx of all fields
  output$download_all <- downloadHandler(
    filename = function() {
      paste("all_fields", ".xlsx", sep = "")
    },
    content = function(file) {
      df <- data.frame(tables()$all_fields)
      write.xlsx(df, file, row.names = FALSE)
    },
    contentType = "application/vnd.ms-excel"
  )
  
  # Downloadable xslx of gene fields
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
  
  # Downloadable xslx of virulence hits
  output$download_vir_hits <- downloadHandler(
    filename = function() {
      paste("virulence_hits", ".xlsx", sep = "")
    },
    content = function(file) {
      df <- data.frame(tables()$virulence_hits)
      write.xlsx(df, file, row.names = FALSE)
    },
    contentType = "application/vnd.ms-excel"
  )
  
  # Downloadable xslx of statistical analysis (chi-square values)
  output$download_stats_ana <- downloadHandler(
    filename = function() {
      paste("stats_analysis", ".xlsx", sep = "")
    },
    content = function(file) {
      df <- data.frame(tables()$stats_analysis)
      write.xlsx(df, file, row.names = FALSE)
    },
    contentType = "application/vnd.ms-excel"
  )
  
  # Call renderPlot for each one. Plots are only actually generated when they
  # are visible on the web page.
  observe({
    displayPlots(tables()$plots$spec_plots, "spec_plot", output)
    displayPlots(tables()$plots$sero_plots, "sero_plot", output)
    displayPlots(tables()$plots$seq_type_plots, "seq_type_plot", output)
  })
})
