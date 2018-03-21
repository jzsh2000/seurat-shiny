#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
# Wed Jan 17 21:49:01 2018 ------------------------------


library(shiny)
library(shinyjs)
library(tidyverse)

if (file.exists('config.txt')) {
    source('config.txt')
}
resource_list <- read_csv('data/resource_list.csv',
                          col_types = 'ccccd')

# Define UI for application that draws a histogram
shinyUI(fluidPage(

  tags$head(tags$link(rel = "shortcut icon", href = "dc.ico")),

  # Application title
  titlePanel(ifelse(exists('app_title'), app_title, "single-cell RNA-seq data visualization")),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
        selectizeInput(inputId = 'dataset',
                       label = 'Dataset',
                       choices = c('(none)' = 'none',
                                   deframe(resource_list[,c('description',
                                                            'label')]))),
        checkboxInput(inputId = 'cb_cellranger',
                      label = 'Use original t-SNE',
                      value = FALSE),
        conditionalPanel(
            'input.tabset_main == "gene expression"',
            # checkboxInput(inputId = 'cb_label',
            #               label = 'Label cluster center',
            #               value = TRUE),
            checkboxInput(inputId = 'cb_showsize',
                          label = 'Show cluster size',
                          value = FALSE)
        ),
        disabled(
            sliderInput(inputId = 'resolution',
                        label = 'Cluster resolution',
                        min = 0.1, max = 1.5, value = 0.8,
                        step = 0.1)
        ),
        conditionalPanel(
            'input.tabset_main == "gene expression"',
            # checkboxInput(inputId = 'cb_label',
            #               label = 'Label cluster center',
            #               value = TRUE),
            textInput(inputId = 'tx_gene',
                      label = 'Gene name',
                      placeholder = 'Your awesome gene')
        ),
        conditionalPanel(
            'input.tabset_main == "co-expression"',
            selectizeInput(inputId = 'cluster_id',
                           label = 'Use cluster',
                           choices = NULL,
                           selected = NULL,
                           multiple = TRUE
                           ),
            textInput(inputId = 'tx_gene1',
                      label = 'Gene name',
                      placeholder = 'Your awesome gene'),
            textInput(inputId = 'tx_gene2',
                      label = 'Gene name',
                      placeholder = 'Your awesome gene')
        ),
        conditionalPanel(
            'input.tabset_main == "signature"',
            selectizeInput(inputId = 'sig_cluster_1',
                           label = 'Use cluster',
                           choices = '(none)',
                           multiple = FALSE),
            selectizeInput(inputId = 'sig_cluster_2',
                           label = 'Compare to',
                           choices = c('all other cells'),
                           selected = NULL,
                           multiple = FALSE)
        )
    ),

    # Show a plot of the generated distribution
    mainPanel(
        tabsetPanel(id = 'tabset_main', type = 'pills',
                    tabPanel(
                        title = 'gene expression',
                        plotOutput('plot_gene_expr')
                    ),
                    tabPanel(
                        title = 'co-expression',
                        verbatimTextOutput('coefficient'),
                        plotOutput('plot_gene_expr2')
                    ),
                    tabPanel(
                        title = 'signature',
                        DT::dataTableOutput('table_sig_gene')
                    )
        )

    )
  )
))
