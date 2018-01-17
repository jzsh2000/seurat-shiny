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
library(here)

resource_list <- read_csv('data/resource_list.csv',
                          col_types = 'cccc')

# Define UI for application that draws a histogram
shinyUI(fluidPage(

  # Application title
  titlePanel("single-cell RNA-seq data visualization"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
        selectizeInput(inputId = 'dataset',
                       label = 'Dataset',
                       choices = c('(none)' = 'none',
                                   deframe(resource_list[,c('description',
                                                            'label')]))),
        conditionalPanel(
            'input.tabset_main == "gene expression"',
            checkboxInput(inputId = 'cb_label',
                          label = 'Label cluster center',
                          value = TRUE),
            checkboxInput(inputId = 'cb_cellranger',
                          label = 'Use cellranger t-SNE',
                          value = FALSE),
            disabled(
                selectizeInput(inputId = 'resolution',
                               label = 'Cluster resolution',
                               choices = c())
            ),
            textInput(inputId = 'tx_gene',
                      label = 'Gene name',
                      placeholder = 'Your awesome gene')
        )
    ),

    # Show a plot of the generated distribution
    mainPanel(
        tabsetPanel(id = 'tabset_main', type = 'pills',
                    tabPanel(
                        title = 'gene expression',
                        plotOutput('plot_gene_expr')
                    )
                    )

    )
  )
))
