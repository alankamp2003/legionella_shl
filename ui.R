library(shiny)
library(shinyjs)
library(DT)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  useShinyjs(),
  # Application title
  titlePanel("Virulence Data Analysis"),
  
  # Sidebar with several inputs for uploading and filtering data
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "Choose virulence data file",
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
      numericInput("scatter_noise", "Scatter Noise (+/-)",
                   min = 0.0, max = 100,
                   value = 0.0, step = 0.1),
      tags$h4("Virulence hits breakdown by"),
      fluidRow(column(4, checkboxInput("vir_spec", label = "Species", value = TRUE)),
               column(4, checkboxInput("vir_sero", label = "Serogroup", value = FALSE)),
               column(4, checkboxInput("vir_seq_type", label = "Seq. Type", value = FALSE))),
      actionButton(inputId = "gen_output",
                   label = "Generate output")
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
                 DT::dataTableOutput('gene_fields')),
        tabPanel("Virulence Hits",
                 tags$p(),
                 downloadButton('download_vir_hits'),
                 tags$p(),
                 DT::dataTableOutput('virulence_hits')),
        tabPanel("Statistical Analysis",
                 tags$p(),
                 downloadButton('download_stats_ana'),
                 tags$p(),
                 DT::dataTableOutput('stats_analysis')),
        tabPanel("Visualization",
                 #tags$p(),
                 #downloadButton('download_vis'),
                 #tags$p(),
                 tabsetPanel(tabPanel("Gene scatter",
                                      tags$h4("By species and facility", style = "color:purple"),
                                      plotOutput("spec_plot", height = "210px"),
                                      tags$h4("By serogroup and facility", style = "color:purple"),
                                      plotOutput("sero_plot", height = "210px"),
                                      tags$h4("By sequence type and facility", style = "color:purple"),
                                      plotOutput("seq_type_plot", height = "210px")),
                             tabPanel("Gene percentage",
                                      tags$h4("By species", style = "color:purple"), 
                                      uiOutput('spec_charts'),
                                      tags$h4("By serogroup", style = "color:purple"), 
                                      uiOutput('sero_charts'),
                                      tags$h4("By sequence type", style = "color:purple"),
                                      uiOutput('seq_type_charts'))))
      )
    )
  )
))
