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
library(shinythemes)
library(shinyjs)
library(tidyverse)
library(glue)

resource_list <- read_csv('data/resource_list.csv',
                          col_types = 'ccccd')

# Read configure parameters
if (file.exists('config.txt')) {
    source('config.txt')
}
if (!exists('app_title')) app_title = "single-cell RNA-seq data visualization"
if (!exists('panel_name_1')) panel_name_1 = 'Gene expression value'
if (!exists('panel_name_2')) panel_name_2 = 'Cluster signature genes'
if (!exists('panel_name_3')) panel_name_3 = 'Gene co-expression'
if (file.exists('www/About.Rmd')) {
    about_page_path = 'www/About.Rmd'
} else {
    about_page_path = 'www/About.Rmd.template'
}

# Define UI
shinyUI(navbarPage(

  theme = shinytheme('cerulean'),
  title = app_title,
  position = 'fixed-top',
  id = 'dc_sc',
  collapsible = TRUE,
  header = tags$head(tags$link(rel = "shortcut icon",
                               href = "dc.ico"),
                     tags$link(rel = "stylesheet",
                               type = "text/css",
                               href = "style.css")),
  selected = 'datasets',

  tabPanel(
      tags$b('About'), value = 'about',
      column(1),
      column(
          width = 10,
          wellPanel(includeMarkdown(about_page_path))
      ),
      column(1)
  ),

  tabPanel(
      tags$b('Single-cell Datasets'),
      value = 'datasets',
      useShinyjs(),

      sidebarLayout(
          sidebarPanel(
              selectizeInput(inputId = 'dataset',
                             label = 'Dataset [*]',
                             choices = c(' ' = 'none',
                                         deframe(resource_list[,c('description',
                                                                  'label')]))),
              hidden(
                  tags$div(
                      id = 'dat_config',
                      sliderInput(inputId = 'resolution',
                                  label = 'Cluster resolution',
                                  min = 0.1, max = 1.5, value = 0.8,
                                  step = 0.1),
                      hr(),
                      conditionalPanel(
                          glue('input.tabset_main == "{panel_name_1}"'),
                          # checkboxInput(inputId = 'cb_label',
                          #               label = 'Label cluster center',
                          #               value = TRUE),
                          #
                          checkboxInput(inputId = 'cb_cellranger',
                                        label = 'Use original t-SNE',
                                        value = FALSE),
                          checkboxInput(inputId = 'cb_allpt',
                                        label = 'Show all cells',
                                        value = FALSE),
                          checkboxInput(inputId = 'cb_showsize',
                                        label = 'Show cluster size',
                                        value = FALSE),
                          checkboxInput(inputId = 'cb_subset',
                                        label = 'Use subset cluster',
                                        value = FALSE),
                          textInput(inputId = 'tx_gene',
                                    label = 'Gene name',
                                    placeholder = 'Your awesome gene'),
                          hr(),
                          selectizeInput(inputId = 'cluster_id_subset',
                                         label = 'Use cluster',
                                         choices = NULL,
                                         selected = NULL,
                                         multiple = TRUE
                          ),
                          sliderInput(inputId = 'resolution_subset',
                                      label = 'Cluster resolution',
                                      min = 0.1, max = 1.5, value = 0.8,
                                      step = 0.1)
                      ),
                      conditionalPanel(
                          glue('input.tabset_main == "{panel_name_2}"'),
                          selectizeInput(inputId = 'cluster_id',
                                         label = 'Use cluster',
                                         choices = NULL,
                                         selected = NULL,
                                         multiple = TRUE
                          ),
                          textInput(inputId = 'tx_gene1',
                                    label = 'Gene name #1 [*]',
                                    placeholder = 'Your awesome gene'),
                          textInput(inputId = 'tx_gene2',
                                    label = 'Gene name #2 [*]',
                                    placeholder = 'Your awesome gene')
                      ),
                      conditionalPanel(
                          glue('input.tabset_main == "{panel_name_3}"'),
                          selectizeInput(inputId = 'sig_cluster_1',
                                         label = 'Use cluster [*]',
                                         choices = '',
                                         multiple = FALSE),
                          selectizeInput(inputId = 'sig_cluster_2',
                                         label = 'Compare to',
                                         choices = '(All other cells)',
                                         # selected = NULL,
                                         multiple = FALSE),
                          radioButtons(inputId = 'marker_pos',
                                       label = 'Find markers',
                                       inline = TRUE,
                                       choices = c('positive' = 'pos',
                                                   'negative' = 'neg'))
                      )
                  )

              )
          ),

          # Show a plot of the generated distribution
          mainPanel(
              hidden(
                  tags$div(
                      id = 'dat_panel',
                      tabsetPanel(id = 'tabset_main', type = 'tabs',
                                  tabPanel(
                                      title = panel_name_1,
                                      plotOutput('plot_gene_expr')
                                  ),
                                  tabPanel(
                                      title = panel_name_2,
                                      tags$p(id = 'warning_info'),
                                      DT::dataTableOutput('table_sig_gene')
                                  ),
                                  tabPanel(
                                      title = panel_name_3,
                                      # verbatimTextOutput('coefficient'),
                                      plotOutput('plot_gene_expr2')
                                  )
                      )
                  )
              )
          )
      )
  )
))
