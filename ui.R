library(shiny)
library(shinyjs)
library(DT)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  useShinyjs(),
  # Application title
  titlePanel("Legionella Data Analysis"),
  
  # Sidebar with several inputs for uploading and filtering data
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "Choose legionella data file",
                accept = c(
                  "application/vnd.ms-excel",
                  "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                  ".xls",
                  ".xlsx")
      ),
      tags$h4("Filter by"),
      selectInput("species", "Species", multiple = TRUE,
                  c("All" = "")),
      fluidRow(column(6, selectInput("serogroup", "Serogroups", multiple = TRUE,
                                     c("All" = ""))),
               column(6, selectInput("seq_type", "Seq. types", multiple = TRUE,
                                     c("All" = "")))),
      selectInput("gene", "Genes", multiple = TRUE,
                  c("All" = "")),
      fluidRow(column(6, numericInput("coverage", "Seq Total % Coverage",
                  min = 0.0, max = 100,
                  value = 0.0, step = 0.5)),
               column(6, numericInput("identity", "Seq Weighted % Identity",
                  min = 0.0, max = 100,
                  value = 0.0, step = 0.5))),
      actionButton(inputId = "gen_tables",
                   label = "Generate tables")
    ),

    # Show a tables for all fields and gene fields
    mainPanel(
      tabsetPanel(
        id = 'dataset',
        tabPanel("All fields", 
                 tags$p(), 
                 downloadButton('download_all'),
                 tags$p(),
                 DT::dataTableOutput('all_fields')),
        tabPanel("Gene fields",
                 tags$p(),
                 downloadButton('download_genes'),
                 tags$p(),
                 DT::dataTableOutput('gene_fields'))
      )
    )
  )
))
