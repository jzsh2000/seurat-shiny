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
                          col_types = 'ccccc')

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
    get_dataset <- reactive({
        if (input$dataset == 'none') {
            list(species = NULL,
                 rdat = NULL,
                 rdat_tsne = NULL
                 )
        } else {
            resource = resource_list %>%
                filter(label == input$dataset)
            rdat = read_rds(file.path('data/robj', resource$robj))
            rdat_tsne = read_rds(file.path('data/robj', resource$tsne))

            resolution_list = str_subset(colnames(rdat@meta.data), '^res\\.') %>%
                str_extract('(?<=res.).*') %>%
                as.numeric() %>%
                sort() %>%
                as.character()

            updateSelectizeInput(session, 'resolution',
                                 choices = resolution_list)
            enable('resolution')

            list(species = resource$species,
                 rdat = rdat,
                 rdat_tsne = rdat_tsne
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
        res_name = paste0('res.', input$resolution)
        if (!is.null(get_dataset()$rdat)) {
            if (get_input_gene() != '' &&
                (get_input_gene() %in% rownames(get_dataset()$rdat@data))) {

                if (input$cb_cellranger) {
                    cluster_dat <- get_dataset()$rdat_tsne %>%
                        inner_join(get_dataset()$rdat@meta.data %>%
                                       rownames_to_column('Barcode') %>%
                                       as_data_frame() %>%
                                       dplyr::select(Barcode, one_of(res_name)),
                                   by = 'Barcode')

                    plot_1 <- ggplot(data = cluster_dat,
                           mapping = aes_string(x = 'tSNE_1', y = 'tSNE_2',
                                                color = res_name)) +
                        geom_point(size = 1) +
                        geom_text(cluster_dat %>%
                                      group_by_(res_name) %>%
                                      summarise(tSNE_1 = mean(tSNE_1),
                                                tSNE_2 = mean(tSNE_2)),
                                  mapping = aes_string(x = 'tSNE_1', y = 'tSNE_2', label = res_name),
                                  color = 'black') +
                        coord_fixed() +
                        theme_bw() +
                        theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              legend.position = "none")

                    plot_2 <- ggplot(cluster_dat %>%
                                         inner_join(
                                             enframe(get_dataset()$rdat@data[get_input_gene(),], name = 'Barcode', value = get_input_gene()),
                                             by = 'Barcode'
                                         ),
                                     aes_string(
                                         x = 'tSNE_1',
                                         y = 'tSNE_2',
                                         color = get_input_gene()
                                         )
                                     ) +
                        scale_colour_gradient(low = 'grey', high = 'blue') +
                        geom_point(size = 1) +
                        coord_fixed() +
                        theme_bw() +
                        theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              legend.position = "none")
                } else {
                    plot_1 <- Seurat::TSNEPlot(get_dataset()$rdat,
                                               no.legend = TRUE,
                                               do.label = TRUE,
                                               group.by = res_name) +
                        coord_fixed()

                    plot_2 <- Seurat::FeaturePlot(get_dataset()$rdat,
                                                  features.plot = get_input_gene(),
                                                  cols.use = c('grey', 'blue'),
                                                  do.return = TRUE)[[1]] +
                        theme_bw() +
                        theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              legend.position = "none") +
                        coord_fixed()
                }

                plot_3 <- Seurat::VlnPlot(get_dataset()$rdat,
                                          group.by = res_name,
                                          features.plot = get_input_gene())

                plot_grid(
                    plot_grid(plot_1, plot_2, align = 'h'),
                    plot_3, ncol = 1,
                    rel_heights = c(3, 2)
                )
            } else {
                if (input$cb_cellranger) {
                    cluster_dat <- get_dataset()$rdat_tsne %>%
                        inner_join(get_dataset()$rdat@meta.data %>%
                                       rownames_to_column('Barcode') %>%
                                       as_data_frame() %>%
                                       dplyr::select(Barcode, one_of(res_name)),
                                   by = 'Barcode')

                    ggplot(data = cluster_dat,
                           mapping = aes_string(x = 'tSNE_1', y = 'tSNE_2',
                                                color = res_name)) +
                        geom_point(size = 1) +
                        geom_text(cluster_dat %>%
                                      group_by_(res_name) %>%
                                      summarise(tSNE_1 = mean(tSNE_1),
                                                tSNE_2 = mean(tSNE_2)),
                                  mapping = aes_string(x = 'tSNE_1', y = 'tSNE_2', label = res_name),
                                  color = 'black',
                                  size = 10) +
                        coord_fixed() +
                        theme_bw() +
                        theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              legend.position = "none")
                } else {
                    Seurat::TSNEPlot(get_dataset()$rdat,
                                     no.legend = TRUE,
                                     do.label = TRUE,
                                     label.size = 10,
                                     group.by = paste0('res.',
                                                       input$resolution)
                    ) +
                        coord_fixed()
                }
            }
        }
    # }, width = 800, height = 600)
    }, height = function() {
        session$clientData$output_plot_gene_expr_width
    })
})
