#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyjs)
library(tidyverse)
library(Seurat)
library(DT)
library(cowplot)

resource_list <- read_csv('data/resource_list.csv',
                          col_types = 'cccc')

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
    get_dataset <- reactive({
        if (input$dataset == 'none') {
            list(species = NULL,
                 rdat = NULL
                 )
        } else {
            resource = resource_list %>%
                filter(label == input$dataset)
            rdat = read_rds(file.path('data/robj', resource$robj))
            rdat_tsne = read_rds(file.path('data/robj', resource$projection))

            resolution_list = str_subset(colnames(rdat@meta.data), '^res\\.') %>%
                str_extract('(?<=res.).*') %>%
                as.numeric() %>%
                sort() %>%
                as.character()

            updateSelectizeInput(session, 'resolution',
                                 choices = resolution_list)
            enable('resolution')

            list(species = resource$species,
                 rdat = rdat
                 )
        }
    }) %>% debounce(1000)

    get_input_gene <- reactive({
        if (is.null(get_dataset()$species)) {
            ''
        } else {
            if (input$tx_gene %in% rownames(rdat@data)) {
                input$tx_gene
            } else {
                if (get_dataset()$species == 'human') {
                    c(limma::alias2Symbol(stringr::str_to_upper(input$tx_gene),
                                          species = 'Hs'), '')[1]
                } else if (get_dataset()$species == 'mouse') {
                    c(limma::alias2Symbol(stringr::str_to_title(input$tx_gene),
                                          species = 'Mm'), '')[1]
                } else {
                    ''
                }
            }
        }
    }) %>% debounce(1500)

    output$plot_tsne <- renderPlot({
        if (!is.null(get_dataset()$rdat)) {

        }
    })

    output$plot_gene_expr <- renderPlot({
        if (!is.null(get_dataset()$rdat)) {
            if (get_input_gene() != '' &&
                (get_input_gene() %in% rownames(get_dataset()$rdat@data))) {

                plot_1 <- Seurat::TSNEPlot(get_dataset()$rdat,
                                           no.legend = input$cb_label,
                                           do.label = input$cb_label,
                                           group.by = paste0('res.', input$resolution)) +
                    coord_fixed()

                plot_2 <- Seurat::FeaturePlot(get_dataset()$rdat,
                                    features.plot = get_input_gene(),
                                    cols.use = c('grey', 'blue'),
                                    do.return = TRUE)[[1]] +
                    coord_fixed()

                plot_3 <- Seurat::VlnPlot(get_dataset()$rdat,
                                          group.by = paste0('res.', input$resolution),
                                          features.plot = get_input_gene())

                plot_grid(
                    plot_grid(plot_1, plot_2, align = 'h'),
                    plot_3, ncol = 1,
                    rel_heights = c(3, 2)
                )
            } else {
                Seurat::TSNEPlot(get_dataset()$rdat,
                                 no.legend = input$cb_label,
                                 do.label = input$cb_label,
                                 label.size = 10,
                                 group.by = paste0('res.',
                                                  input$resolution)
                                 ) +
                    coord_fixed()
            }
        }
    # }, width = 800, height = 600)
    }, height = function() {
        session$clientData$output_plot_gene_expr_width
    })
})
