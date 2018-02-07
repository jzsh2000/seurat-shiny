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
    plot_height_func <- function(scale_factor) {
        return(function() {
            session$clientData$output_plot_gene_expr_width * scale_factor
        })
    }

    # get standard gene name
    gene_name_std <- function(gene_name, species = 'human') {
        if (species == 'human') {
            input_gene = c(limma::alias2Symbol(stringr::str_to_upper(gene_name),
                                               species = 'Hs'), '')[1]
        } else if (species == 'mouse') {
            input_gene = c(limma::alias2Symbol(stringr::str_to_title(gene_name),
                                               species = 'Mm'), '')[1]
        } else {
            input_gene = ''
        }
        input_gene
    }

    get_dataset <- reactive({
        if (input$dataset == 'none') {
            list(species = NULL,
                 rdat = NULL,
                 rdat_tsne = NULL
                 )
        } else {
            withProgress(message = 'Load seurat object',
                         detail = 'Locate RDS file path',
                         value = 0, {
                             resource = resource_list %>%
                                 filter(label == input$dataset)
                             incProgress(0.1, message = 'Read RDS file')

                             rdat = read_rds(file.path('data/robj', resource$robj))
                             incProgress(0.6, message = 'Read t-SNE coordinates')

                             rdat_tsne = read_rds(file.path('data/robj', resource$tsne))
                             incProgress(0.2, message = 'Get resolution list')

                             resolution_list = str_subset(colnames(rdat@meta.data), '^res\\.') %>%
                                 str_extract('(?<=res.).*') %>%
                                 as.numeric() %>%
                                 sort() %>%
                                 as.character()

                             updateSelectizeInput(session, 'resolution',
                                                  choices = resolution_list)
                             enable('resolution')

                             setProgress(value = 1, message = 'Finish!')
                         })

            list(species = resource$species,
                 rdat = rdat,
                 rdat_tsne = rdat_tsne
                 )
        }
    }) %>% debounce(1000)

    observe(if (!is.null(input$resolution) &&
                !is.null(get_dataset()$rdat)) {
        res_name = paste0('res.', input$resolution)
        updateSelectizeInput(session, 'cluster_id', choices = as.character(sort(as.integer(unique(get_dataset()$rdat@meta.data[[res_name]])))), selected = '0')
    })

    get_input_gene <- reactive({
        input_gene = ''

        if (!is.null(get_dataset()$species)) {
            if (input$tx_gene %in% rownames(get_dataset()$rdat@data)) {
                input_gene = input$tx_gene
            } else {
                input_gene = gene_name_std(input$tx_gene, get_dataset()$species)
            }
        }

        if (input_gene != '' &&
            (input_gene %in% rownames(get_dataset()$rdat@data))) {
            input_gene
        } else {
            ''
        }
    }) %>% debounce(1500)

    get_input_gene1 <- reactive({
        input_gene = ''

        if (!is.null(get_dataset()$species)) {
            if (input$tx_gene1 %in% rownames(get_dataset()$rdat@data)) {
                input_gene = input$tx_gene1
            } else {
                input_gene = gene_name_std(input$tx_gene1, get_dataset()$species)
            }
        }

        if (input_gene != '' &&
            (input_gene %in% rownames(get_dataset()$rdat@data))) {
            input_gene
        } else {
            ''
        }
    }) %>% debounce(1500)

    get_input_gene2 <- reactive({
        input_gene = ''

        if (!is.null(get_dataset()$species)) {
            if (input$tx_gene2 %in% rownames(get_dataset()$rdat@data)) {
                input_gene = input$tx_gene2
            } else {
                input_gene = gene_name_std(input$tx_gene2, get_dataset()$species)
            }
        }

        if (input_gene != '' &&
            (input_gene %in% rownames(get_dataset()$rdat@data))) {
            input_gene
        } else {
            ''
        }
    }) %>% debounce(1500)

    get_cluster_dat_cellranger <- reactive({
        res_name = paste0('res.', input$resolution)
        cluster_dat <- get_dataset()$rdat_tsne %>%
            inner_join(get_dataset()$rdat@meta.data %>%
                           rownames_to_column('Barcode') %>%
                           as_data_frame() %>%
                           dplyr::select(Barcode, one_of(res_name)),
                       by = 'Barcode') %>%
            dplyr::rename_(cluster = res_name)
        cluster_dat
    })

    get_cluster_dat_seurat <- reactive({
        res_name = paste0('res.', input$resolution)
        cluster_dat <- GetDimReduction(get_dataset()$rdat,
                                       reduction.type = 'tsne',
                                       slot = 'cell.embeddings') %>%
            as.data.frame() %>%
            rownames_to_column(var = 'Barcode') %>%
            as_data_frame() %>%
            inner_join(get_dataset()$rdat@meta.data %>%
                           rownames_to_column('Barcode') %>%
                           as_data_frame() %>%
                           dplyr::select(Barcode, one_of(res_name)),
                       by = 'Barcode') %>%
            dplyr::rename_(cluster = res_name)
        cluster_dat
    })

    get_tsne_plot <- reactive({
        res_name = paste0('res.', input$resolution)
        label_column = 'cluster'
        if (input$cb_showsize) {
            label_column = 'cluster_with_size'
        }

        if (input$cb_cellranger) {
            plot_dat = get_cluster_dat_cellranger()
        } else {
            plot_dat = get_cluster_dat_seurat()
        }

        ggplot(data = plot_dat,
               mapping = aes(x = tSNE_1, y = tSNE_2, color = cluster)) +
            geom_point(size = 1) +
            geom_text(plot_dat %>%
                          group_by(cluster) %>%
                          summarise(tSNE_1 = mean(tSNE_1),
                                    tSNE_2 = mean(tSNE_2),
                                    cluster_size = n()) %>%
                          mutate(cluster_with_size = paste0('Cluster ', cluster, ' (', cluster_size, ')')),
                      mapping = aes_string(x = 'tSNE_1', y = 'tSNE_2',
                                           label = label_column),
                      color = 'black',
                      fontface = "bold") +
            coord_fixed() +
            theme_bw() +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  legend.position = "none")
    })

    output$plot_gene_expr <- renderPlot({
        res_name = paste0('res.', input$resolution)
        if (!is.null(get_dataset()$rdat)) {
            if (get_input_gene() != '' &&
                (get_input_gene() %in% rownames(get_dataset()$rdat@data))) {

                plot_1 <- get_tsne_plot()
                if (input$cb_cellranger) {
                    plot_dat <- get_cluster_dat_cellranger()
                } else {
                    plot_dat <- get_cluster_dat_seurat()
                }

                plot_dat  %>%
                    inner_join(
                        enframe(get_dataset()$rdat@data[get_input_gene(),], name = 'Barcode', value = get_input_gene()),
                        by = 'Barcode'
                    )

                plot_2 <- ggplot(plot_dat,
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

                plot_3 <- Seurat::VlnPlot(get_dataset()$rdat,
                                          group.by = res_name,
                                          features.plot = get_input_gene())

                plot_grid(
                    plot_grid(plot_1, plot_2, align = 'h'),
                    plot_3, ncol = 1,
                    rel_heights = c(3, 2)
                )
            } else {
                get_tsne_plot()
            }
        }
    # }, width = 800, height = 600)
    }, height = plot_height_func(0.7))

    shared_data <- reactiveValues(cor = NULL)

    output$plot_gene_expr2 <- renderPlot({
        res_name = paste0('res.', input$resolution)
        if (!is.null(get_dataset()$rdat) &&
            get_input_gene1() != '' &&
            (get_input_gene1() %in% rownames(get_dataset()$rdat@data)) &&
            get_input_gene2() != '' &&
            (get_input_gene2() %in% rownames(get_dataset()$rdat@data)) &&
            get_input_gene1() != get_input_gene2()) {

            valid_cells = rownames(get_dataset()$rdat@meta.data)[get_dataset()$rdat@meta.data[[res_name]] %in% input$cluster_id]

            plot_dat <- inner_join(
                enframe(get_dataset()$rdat@data[get_input_gene1(),], name = 'Barcode', value = get_input_gene1()),
                enframe(get_dataset()$rdat@data[get_input_gene2(),], name = 'Barcode', value = get_input_gene2()),
                by = 'Barcode'
            ) %>%
                filter(Barcode %in% valid_cells)

            shared_data$cor = cor(plot_dat[[get_input_gene1()]],
                                  plot_dat[[get_input_gene2()]])
            ggplot(plot_dat,
                aes_string(
                    x = get_input_gene1(),
                    y = get_input_gene2()
                    )
            ) +
            geom_point(size = 1) +
            coord_fixed() +
            theme_bw() +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  legend.position = "none")

        }
    }, height = plot_height_func(0.6))

    output$coefficient <- renderText({
        paste('Pearson correlation coeffient:', shared_data$cor)
    })
})

