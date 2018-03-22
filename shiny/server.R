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
library(magrittr)
library(Seurat)
library(DT)
library(cowplot)

resource_list <- read_csv('data/resource_list.csv',
                          col_types = 'ccccd')

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

    dataset_info <- reactiveValues(
        species = NULL,
        rdat = NULL,
        rdat_tsne = NULL
    )

    update_resolution <- function(resolution) {
        res_name = paste0('res.', resolution)
        if (!is.null(dataset_info$rdat) &&
            !(res_name %in% colnames(dataset_info$rdat@meta.data))) {
            print('update resolution...')
            dataset_info$rdat = FindClusters(
                object = dataset_info$rdat,
                dims.use = 1:15,
                print.output = FALSE,
                resolution = input$resolution
            )
            print('done...')
        }
    }

    observeEvent(input$dataset, {
        if (input$dataset == 'none') {
            dataset_info$species = NULL
            dataset_info$rdat = NULL
            dataset_info$rdat_tsne = NULL
            get_sig_gene$table = data_frame(gene = character())

            updateSelectizeInput(session,
                                 inputId = 'sig_cluster_1',,
                                 choices = '(none)',
                                 selected = '(none)')
            updateSelectizeInput(session,
                                 inputId = 'sig_cluster_2',,
                                 choices = 'all other cells',
                                 selected = 'all other cells')
        } else {
            withProgress(message = 'Load seurat object',
                         detail = 'Locate RDS file path',
                         value = 0, {
                             resource = resource_list %>%
                                 filter(label == input$dataset)
                             incProgress(0.1, message = 'Read RDS file')

                             # print(file.path('data', resource$data_dir, paste0(resource$data_dir, '.rds')))
                             rdat = read_rds(file.path('data', resource$data_dir, paste0(resource$data_dir, '.rds')))
                             incProgress(0.6, message = 'Read t-SNE coordinates')

                             rdat_tsne = read_csv(file.path('data', resource$data_dir, 'projection.csv'), col_types = 'cdd') %>%
                                 mutate(Barcode = str_extract(Barcode, '^[^-]+')) %>%
                                 rename(tSNE_1 = `TSNE-1`, tSNE_2 = `TSNE-2`)

                             incProgress(0.2, message = 'Get resolution list')

                             # resolution_list = str_subset(colnames(rdat@meta.data), '^res\\.') %>%
                             #     str_extract('(?<=res.).*') %>%
                             #     as.numeric() %>%
                             #     sort() %>%
                             #     as.character()

                             default_resolution = round(resource$default_resolution, digits = 1)
                             if (!is.na(default_resolution) &&
                                 default_resolution >= 0.1 &&
                                 default_resolution <= 1.5) {
                                 updateSliderInput(session, 'resolution',
                                                   value = default_resolution)
                             }

                             enable('resolution')

                             setProgress(value = 1, message = 'Finish!')
                         })

            dataset_info$species = resource$species
            dataset_info$rdat = rdat
            dataset_info$rdat_tsne = rdat_tsne
            get_sig_gene$table = data_frame(gene = character())
        }
    })

    observe(if (!is.null(input$resolution) &&
                !is.null(dataset_info$rdat)) {
        res_name = paste0('res.', input$resolution)
        update_resolution(input$resolution)
        res_choices = as.character(sort(as.integer(unique(dataset_info$rdat@meta.data[[res_name]]))))

        updateSelectizeInput(session, 'cluster_id',
                             choices = res_choices,
                             selected = '0')
        updateSelectizeInput(session, 'sig_cluster_1',
                             choices = c('(none)', res_choices),
                             selected = '(none)')
        updateSelectizeInput(session, 'sig_cluster_2',
                             choices = c('all other cells', res_choices),
                             selected = 'all other cells')
    })

    get_input_gene <- reactive({
        input_gene = ''

        if (!is.null(dataset_info$species)) {
            if (input$tx_gene %in% rownames(dataset_info$rdat@data)) {
                input_gene = input$tx_gene
            } else {
                input_gene = gene_name_std(input$tx_gene, dataset_info$species)
            }
        }

        if (input_gene != '' &&
            (input_gene %in% rownames(dataset_info$rdat@data))) {
            input_gene
        } else {
            ''
        }
    }) %>% debounce(1500)

    get_input_gene1 <- reactive({
        input_gene = ''

        if (!is.null(dataset_info$species)) {
            if (input$tx_gene1 %in% rownames(dataset_info$rdat@data)) {
                input_gene = input$tx_gene1
            } else {
                input_gene = gene_name_std(input$tx_gene1, dataset_info$species)
            }
        }

        if (input_gene != '' &&
            (input_gene %in% rownames(dataset_info$rdat@data))) {
            input_gene
        } else {
            ''
        }
    }) %>% debounce(1500)

    get_input_gene2 <- reactive({
        input_gene = ''

        if (!is.null(dataset_info$species)) {
            if (input$tx_gene2 %in% rownames(dataset_info$rdat@data)) {
                input_gene = input$tx_gene2
            } else {
                input_gene = gene_name_std(input$tx_gene2, dataset_info$species)
            }
        }

        if (input_gene != '' &&
            (input_gene %in% rownames(dataset_info$rdat@data))) {
            input_gene
        } else {
            ''
        }
    }) %>% debounce(1500)

    get_cluster_dat_cellranger <- reactive({
        res_name = paste0('res.', input$resolution)
        cluster_dat <- dataset_info$rdat_tsne %>%
            left_join(dataset_info$rdat@meta.data %>%
                           rownames_to_column('Barcode') %>%
                           as_data_frame() %>%
                           dplyr::select(Barcode, one_of(res_name)),
                       by = 'Barcode') %>%
            dplyr::rename_(cluster = res_name)
        cluster_dat
    })

    get_cluster_dat_seurat <- reactive({
        res_name = paste0('res.', input$resolution)
        cluster_dat <- GetDimReduction(dataset_info$rdat,
                                       reduction.type = 'tsne',
                                       slot = 'cell.embeddings') %>%
            as.data.frame() %>%
            rownames_to_column(var = 'Barcode') %>%
            as_data_frame() %>%
            inner_join(dataset_info$rdat@meta.data %>%
                           rownames_to_column('Barcode') %>%
                           as_data_frame() %>%
                           dplyr::select(Barcode, one_of(res_name)),
                       by = 'Barcode') %>%
            dplyr::rename_(cluster = res_name)
        cluster_dat
    })

    get_tsne_plot <- reactive({
        update_resolution(input$resolution)
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
        update_resolution(input$resolution)

        if (!is.null(dataset_info$rdat)) {
            if (get_input_gene() != '' &&
                (get_input_gene() %in% rownames(dataset_info$rdat@data))) {

                plot_1 <- get_tsne_plot()
                if (input$cb_cellranger) {
                    plot_dat <- get_cluster_dat_cellranger()
                } else {
                    plot_dat <- get_cluster_dat_seurat()
                }

                plot_dat  %<>%
                    inner_join(
                        enframe(dataset_info$rdat@data[get_input_gene(),], name = 'Barcode', value = 'expr'),
                        by = 'Barcode'
                    )

                plot_2 <- ggplot(plot_dat,
                                 aes(
                                     x = tSNE_1,
                                     y = tSNE_2,
                                     color = expr
                                 )
                ) +
                    scale_colour_gradient(low = 'grey', high = 'blue') +
                    geom_point(size = 1) +
                    coord_fixed() +
                    theme_bw() +
                    theme(panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          legend.position = "none")

                plot_3 <- Seurat::VlnPlot(dataset_info$rdat,
                                          group.by = res_name,
                                          features.plot = get_input_gene(),
                                          do.return = TRUE)

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
        if (!is.null(dataset_info$rdat) &&
            get_input_gene1() != '' &&
            (get_input_gene1() %in% rownames(dataset_info$rdat@data)) &&
            get_input_gene2() != '' &&
            (get_input_gene2() %in% rownames(dataset_info$rdat@data)) &&
            get_input_gene1() != get_input_gene2()) {

            valid_cells = rownames(dataset_info$rdat@meta.data)[dataset_info$rdat@meta.data[[res_name]] %in% input$cluster_id]

            plot_dat <- inner_join(
                enframe(dataset_info$rdat@data[get_input_gene1(),], name = 'Barcode', value = get_input_gene1()),
                enframe(dataset_info$rdat@data[get_input_gene2(),], name = 'Barcode', value = get_input_gene2()),
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

    get_sig_gene <- reactiveValues(table = NULL)
    find_marker_gene <- function() {
        withProgress(message = 'Find marker gene',
                     detail = 'collect cell id in cluster',
                     value = 0, {
            res_name = paste0('res.', input$resolution)
            cell_1 = rownames(dataset_info$rdat@meta.data)[dataset_info$rdat@meta.data[[res_name]] %in% input$sig_cluster_1]

            incProgress(amount = 0.2, message = 'run program')
            if (input$sig_cluster_2[1] == 'all other cells') {
                output_df = FindMarkers(dataset_info$rdat,
                                        ident.1 = cell_1,
                                        ident.2 = NULL,
                                        test.use = 'roc',
                                        min.pct = 0.25,
                                        only.pos = TRUE
                                        )
            } else {
                cell_2 = rownames(dataset_info$rdat@meta.data)[dataset_info$rdat@meta.data[[res_name]] %in% input$sig_cluster_2]
                output_df = FindMarkers(dataset_info$rdat,
                                        ident.1 = cell_1,
                                        ident.2 = cell_2,
                                        test.use = 'roc',
                                        min.pct = 0.25,
                                        only.pos = TRUE
                )
            }

            incProgress(amount = 0.7, message = 'create output dataframe')
            output_df %<>%
                rownames_to_column(var = 'gene') %>%
                as_data_frame() %>%
                dplyr::select(-p_val_adj) %>%
                dplyr::filter(myAUC >= 0.7) %>%
                arrange(desc(myAUC))

            setProgress(value = 1, message = 'Finish!')
        })
        output_df
    }

    get_sig_cluster_input <- reactive({
        list(cluster1 = input$sig_cluster_1,
             cluster2 = input$sig_cluster_2)
    }) %>% debounce(2000)

    observeEvent(get_sig_cluster_input(), {
        if (length(input$sig_cluster_1) > 0 &&
            length(input$sig_cluster_2) > 0 &&
            input$sig_cluster_1[1] != '(none)' &&
            length(intersect(input$sig_cluster_1, input$sig_cluster_2)) == 0) {

            get_sig_gene$table = find_marker_gene()
        }
    })

    output$table_sig_gene <- DT::renderDataTable({
        get_sig_gene$table
    }, selection = 'single')

    observe(
        if (!is.null(input$table_sig_gene_row_last_clicked)) {
            gene_name = get_sig_gene$table$gene[
                input$table_sig_gene_row_last_clicked]

            if (!is.na(gene_name)) {
                updateTextInput(session, "tx_gene", value = gene_name)
            }
        }
    )
})

