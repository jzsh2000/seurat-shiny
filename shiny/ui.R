# Sat Mar  9 23:31:20 2019 ------------------------------

library(shiny)
library(shinythemes)
library(shinyjs)
library(tidyverse)
library(magrittr)
library(glue)
library(jsonlite)

resource.list <- read_json('data/resource_list.json')
source('config.txt')

panel_name_1 = 'Gene expression value'
panel_name_2 = 'Cluster signature genes'
panel_name_3 = 'Gene co-expression'
panel_name_4 = 'Dataset quality'

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
        title = tags$b('About'),
        value = 'about',
        column(1),
        column(
            width = 10,
            wellPanel(includeMarkdown('www/About.Rmd'))
        ),
        column(1)
    ),

    tabPanel(
        title = tags$b('Single-cell Datasets'),
        value = 'datasets',
        useShinyjs(),

        sidebarLayout(
            sidebarPanel(
                selectizeInput(
                    inputId = 'dataset',
                    label = 'Dataset [*]',
                    choices = c(
                        set_names(
                            map_chr(resource.list, ~.$label),
                            map_chr(resource.list, ~.$description)
                        )
                    ),
                    multiple = TRUE,
                    options = list(maxItems = 1)
                ),
                hidden(
                    tags$div(
                        id = 'dat_config',
                        sliderInput(
                            inputId = 'resolution',
                            label = 'Cluster resolution',
                            min = 0, max = 2,
                            value = res_default,
                            step = 0.1
                        ),
                        conditionalPanel(
                            glue('input.tabset_main == "{panel_name_1}"'),
                            selectizeInput(
                                inputId = 'dr_method',
                                label = 'Dimensional Reduction Method',
                                choices = NULL,
                                selected = NULL
                            ),
                            textInput(
                                inputId = 'tx_gene',
                                label = 'Gene name',
                                value = '',
                                placeholder = 'Your awesome gene'
                            ),
                            hr(),
                            selectizeInput(
                                inputId = 'cb_subset',
                                label = 'Use subset cluster',
                                choices = c(
                                    '(None)' = 'none',
                                    'Custom' = 'custom'
                                ),
                                selected = '(None)',
                                multiple = FALSE
                            ),
                            checkboxInput(inputId = 'cb_allpt',
                                          label = 'Show all cells',
                                          value = TRUE),
                            selectizeInput(inputId = 'cluster_id_subset',
                                           label = 'Choose cluster ID',
                                           choices = NULL,
                                           selected = NULL,
                                           multiple = TRUE
                            ),
                            sliderInput(inputId = 'resolution_subset',
                                        label = 'Cluster resolution',
                                        min = 0.1, max = 2.0,
                                        value = res_default,
                                        step = 0.1)
                        ),
                        conditionalPanel(
                            glue('input.tabset_main == "{panel_name_2}"'),
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
                        ),
                        conditionalPanel(
                            glue('input.tabset_main == "{panel_name_3}"'),
                            selectizeInput(inputId = 'cluster_id',
                                           label = 'Use cluster [*]',
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
                        hr(),
                        textOutput(outputId = 'dat_info_text'),
                        conditionalPanel(
                            glue('input.tabset_main == "{panel_name_1}"'),
                            downloadButton(
                                'd_img',
                                label = 'Download PDF image'
                            )
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
                                    ),
                                    tabPanel(
                                        title = panel_name_4,
                                        plotOutput('plot_data_quality')
                                    )
                        )
                    )
                )
            )
        )
    )
))
