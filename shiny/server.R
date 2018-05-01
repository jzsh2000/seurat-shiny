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
library(glue)

## ---------- key configurations

# color for gene name validation: grey (empty) -- red (invalid) -- green (valid)
color_primary = '#CCCCCC'
color_error = '#FF0000'
color_success = '#00FF00'

# default resolution
res_default = 0.8

# plot view scale factor
scale_factor_default = 0.7

## ---------- load

# data.frame for all available dataset
resource_list <- read_csv('data/resource_list.csv',
                          col_types = 'ccccd') %>%
    replace_na(list(default_resolution = res_default))

# gene official symbols and their synonyms
# see: ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/
gene_human <- read_csv('gene/gene-human.csv', col_types = 'cc')
gene_mouse <- read_csv('gene/gene-mouse.csv', col_types = 'cc')

## ----------

# data.frame for signature gene
empty_sig_df <- data_frame(
    gene = character(),
    myAUC = numeric(),
    avg_diff = numeric(),
    power = numeric(),
    avg_logFC = numeric(),
    pct.1 = numeric(),
    pct.2 = numeric()
)

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
    plot_height_func <- function(scale_factor) {
        return(function() {
            session$clientData$output_plot_gene_expr_width * scale_factor
        })
    }

    # get standard gene name
    gene_name_std <- function(gene_name, species = 'human', valid_names = NULL) {
        if (!is.null(valid_names) && (gene_name %in% valid_names)) {
            input_gene = gene_name
        } else {
            if (species == 'human') {
                input_gene = c(limma::alias2Symbol(stringr::str_to_upper(gene_name), species = 'Hs'), '')[1]
            } else if (species == 'mouse') {
                input_gene = c(limma::alias2Symbol(stringr::str_to_title(gene_name), species = 'Mm'), '')[1]
            } else {
                input_gene = ''
            }
        }

        if (!is.null(valid_names) && (!(input_gene %in% valid_names))) {
            input_gene = ''
        }
        input_gene
    }

    # create title for violin plot (show gene symbol and all synonyms)
    gene_name_title <- function(gene_name, species = 'human') {
        shorten_alias <- function(string) {
            alias_full = str_split(string, pattern = '\\|')[[1]]
            if (length(alias_full) > 5) {
                return(paste(alias_full[1:5], collapse = '|'))
            } else {
                return(string)
            }
        }
        if (species == 'human') {
            if (gene_name %in% gene_human$gene) {
                gene_idx = which(gene_name == gene_human$gene)[1]
                gene_alias = shorten_alias(gene_human$alias[gene_idx])
                glue::glue('{gene_name} ({gene_alias})')
            } else {
                gene_name
            }
        } else if (species == 'mouse') {
            if (gene_name %in% gene_mouse$gene) {
                gene_idx = which(gene_name == gene_mouse$gene)[1]
                gene_alias = shorten_alias(gene_mouse$alias[gene_idx])
                # print(gene_alias)
                glue::glue('{gene_name} ({gene_alias})')
            } else {
                gene_name
            }
        } else {
            gene_name
        }
    }

    dataset_info <- reactiveValues(
        species = NULL,
        rdat = NULL,
        # cellranger t-SNE coordinates
        rdat_tsne_cr = NULL,
        # seurat t-SNE coordinates
        rdat_tsne_sr_full = NULL,
        rdat_tsne_sr = NULL,
        resolution = NULL
    )

    get_sig_gene <- reactiveValues(table = NULL)

    observeEvent(input$resolution, {
        req(dataset_info$resolution)
        if (input$resolution >= 0.1 && input$resolution <= 1.5) {
            dataset_info$resolution = input$resolution
        }
    })

    update_resolution <- function(resolution) {
        res_name = glue('res.{resolution}')
        if (!is.null(dataset_info$rdat) &&
            !(res_name %in% colnames(dataset_info$rdat@meta.data))) {
            print('update resolution...')
            withProgress(
                message = 'Find clusters using new reolution',
                detail = 'This may take a while...',
                value = 0, {
                    dataset_info$rdat = FindClusters(
                        object = dataset_info$rdat,
                        dims.use = 1:15,
                        print.output = FALSE,
                        resolution = resolution
                    )
                    setProgress(value = 1)
                }
            )
        }
    }

    observeEvent(input$dataset, {
        if (input$dataset == 'none') {
            shinyjs::hide(id = 'dat_config')
            shinyjs::hide(id = 'dat_panel')

            dataset_info$species = NULL
            dataset_info$rdat = NULL
            dataset_info$rdat_tsne_cr = NULL
            dataset_info$rdat_tsne_sr_full = NULL
            dataset_info$rdat_tsne_sr = NULL
            dataset_info$resolution = NULL
            get_sig_gene$table = empty_sig_df

            updateSelectizeInput(session,
                                 inputId = 'sig_cluster_1',,
                                 choices = NULL,
                                 selected = NULL)
            updateSelectizeInput(session,
                                 inputId = 'sig_cluster_2',,
                                 choices = NULL,
                                 selected = NULL)
        } else {
            shinyjs::show(id = 'dat_config')
            shinyjs::show(id = 'dat_panel')

            withProgress(message = 'Load seurat object',
                         detail = 'Locate RDS file path',
                         value = 0, {
                             resource = resource_list %>%
                                 filter(label == input$dataset)
                             incProgress(0.1, message = 'Read RDS file')

                             rdat = read_rds(file.path('data', resource$data_dir, glue('{resource$data_dir}.rds')))
                             incProgress(0.6, message = 'Read t-SNE coordinates')

                             rdat_tsne_cr = read_csv(file.path('data', resource$data_dir, 'projection.csv'), col_types = 'cdd') %>%
                                 mutate(Barcode = str_extract(Barcode, '^[^-]+')) %>%
                                 rename(tSNE_1 = `TSNE-1`, tSNE_2 = `TSNE-2`)

                             projection_rds_file = file.path('data', resource$data_dir, 'projection.rds')

                             rdat_tsne_sr = GetDimReduction(rdat, reduction.type = 'tsne', slot = 'cell.embeddings') %>%
                                 as.data.frame() %>%
                                 rownames_to_column(var = 'Barcode') %>%
                                 as_data_frame()
                             if (file.exists(projection_rds_file)) {
                                 rdat_tsne_sr_full = read_rds(projection_rds_file)
                             } else {
                                 rdat_tsne_sr_full = rdat_tsne_sr
                             }

                             incProgress(0.2, message = 'Get resolution list')

                             # resolution_list = str_subset(colnames(rdat@meta.data), '^res\\.') %>%
                             #     str_extract('(?<=res.).*') %>%
                             #     as.numeric() %>%
                             #     sort() %>%
                             #     as.character()

                             default_resolution = round(resource$default_resolution, digits = 1)
                             if (!is.na(default_resolution)) {
                                 if (default_resolution < 0.1 || default_resolution > 1.5) {
                                     default_resolution = 0.8
                                 } else {
                                     updateSliderInput(session, 'resolution',
                                                       value = default_resolution)
                                 }
                             } else {
                                 default_resolution = 0.8
                             }

                             setProgress(value = 1, message = 'Finish!')
                         })

            dataset_info$species = resource$species
            dataset_info$rdat = rdat
            dataset_info$rdat_tsne_cr = rdat_tsne_cr
            dataset_info$rdat_tsne_sr_full = rdat_tsne_sr_full
            dataset_info$rdat_tsne_sr = rdat_tsne_sr
            dataset_info$resolution = default_resolution
            get_sig_gene$table = empty_sig_df
            runjs("document.getElementById('warning_info').innerHTML = ''")
        }
    })

    observe(if (!is.null(dataset_info$resolution) &&
                !is.null(dataset_info$rdat)) {
        res_name = glue('res.{dataset_info$resolution}')
        update_resolution(dataset_info$resolution)
        res_choices = as.character(sort(as.integer(unique(dataset_info$rdat@meta.data[[res_name]]))))

        # multiple selection
        updateSelectizeInput(session, 'cluster_id',
                             choices = res_choices)
        # single selection
        updateSelectizeInput(session, 'sig_cluster_1',
                             choices = c('', res_choices),
                             selected = '')
        # single selection
        updateSelectizeInput(session, 'sig_cluster_2',
                             choices = c('(All other cells)', res_choices),
                             selected = '(All other cells)')
    })

    get_input_gene <- reactive({
        input_gene = ''

        if (input$tx_gene == '') {
            runjs(glue("document.getElementById('tx_gene').style.borderColor='{color_primary}'"))
        }
        else if (!is.null(dataset_info$rdat)) {

            if (input$tx_gene %in% rownames(dataset_info$rdat@data)) {
                input_gene = input$tx_gene
                runjs(glue("document.getElementById('tx_gene').style.borderColor='{color_success}'"))
            } else {
                input_gene = gene_name_std(input$tx_gene,
                                           dataset_info$species,
                                           rownames(dataset_info$rdat@data))
                if (input_gene == '') {
                    runjs(glue("document.getElementById('tx_gene').style.borderColor='{color_error}'"))
                } else {
                    runjs(glue("document.getElementById('tx_gene').style.borderColor='{color_success}'"))
                }
            }
        }
        input_gene
    }) %>% debounce(1500)

    get_input_gene1 <- reactive({
        input_gene = ''

        if (input$tx_gene1 == '') {
            runjs(glue("document.getElementById('tx_gene1').style.borderColor='{color_primary}'"))
        }
        else if (!is.null(dataset_info$rdat)) {

            if (input$tx_gene %in% rownames(dataset_info$rdat@data)) {
                input_gene = input$tx_gene1
                runjs(glue("document.getElementById('tx_gene1').style.borderColor='{color_success}'"))
            } else {
                input_gene = gene_name_std(input$tx_gene1,
                                           dataset_info$species,
                                           rownames(dataset_info$rdat@data))
                if (input_gene == '') {
                    runjs(glue("document.getElementById('tx_gene1').style.borderColor='{color_error}'"))
                } else {
                    runjs(glue("document.getElementById('tx_gene1').style.borderColor='{color_success}'"))
                }
            }
        }
        input_gene
    }) %>% debounce(1500)

    get_input_gene2 <- reactive({
        input_gene = ''

        if (input$tx_gene2 == '') {
            runjs(glue("document.getElementById('tx_gene2').style.borderColor='{color_primary}'"))
        }
        else if (!is.null(dataset_info$rdat)) {

            if (input$tx_gene %in% rownames(dataset_info$rdat@data)) {
                input_gene = input$tx_gene2
                runjs(glue("document.getElementById('tx_gene2').style.borderColor='{color_success}'"))
            } else {
                input_gene = gene_name_std(input$tx_gene2,
                                           dataset_info$species,
                                           rownames(dataset_info$rdat@data))
                if (input_gene == '') {
                    runjs(glue("document.getElementById('tx_gene2').style.borderColor='{color_error}'"))
                } else {
                    runjs(glue("document.getElementById('tx_gene2').style.borderColor='{color_success}'"))
                }
            }
        }
        input_gene
    }) %>% debounce(1500)

    get_cluster_dat_cellranger <- reactive({
        res_name = glue('res.{dataset_info$resolution}')
        cluster_dat <- dataset_info$rdat_tsne_cr %>%
            left_join(dataset_info$rdat@meta.data %>%
                           rownames_to_column('Barcode') %>%
                           as_data_frame() %>%
                           dplyr::select(Barcode, one_of(res_name)),
                       by = 'Barcode') %>%
            dplyr::rename_(cluster = res_name)

        if (!input$cb_allpt) {
            cluster_dat %<>%
                filter(!is.na(cluster))
        }
        cluster_dat
    })

    get_cluster_dat_seurat <- reactive({
        res_name = glue('res.{dataset_info$resolution}')
        if (!input$cb_allpt) {
            cluster_dat <- dataset_info$rdat_tsne_sr %>%
                left_join(dataset_info$rdat@meta.data %>%
                              rownames_to_column('Barcode') %>%
                              as_data_frame() %>%
                              dplyr::select(Barcode, one_of(res_name)),
                          by = 'Barcode') %>%
                dplyr::rename_(cluster = res_name)
        } else {
            cluster_dat <- dataset_info$rdat_tsne_sr_full %>%
                left_join(dataset_info$rdat@meta.data %>%
                              rownames_to_column('Barcode') %>%
                              as_data_frame() %>%
                              dplyr::select(Barcode, one_of(res_name)),
                          by = 'Barcode') %>%
                dplyr::rename_(cluster = res_name)
        }
        cluster_dat
    })

    get_tsne_plot <- reactive({
        update_resolution(dataset_info$resolution)
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
            scale_color_discrete(na.value = 'lightgrey') +
            geom_text(plot_dat %>%
                          group_by(cluster) %>%
                          summarise(tSNE_1 = mean(tSNE_1),
                                    tSNE_2 = mean(tSNE_2),
                                    cluster_size = n()) %>%
                          mutate(cluster_with_size = glue('Cluster {cluster} ({cluster_size})')),
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
        res_name = glue('res.{dataset_info$resolution}')
        update_resolution(dataset_info$resolution)

        if (!is.null(dataset_info$rdat)) {
            if (get_input_gene() != '' &&
                (get_input_gene() %in% rownames(dataset_info$rdat@data))) {

                plot_1 <- get_tsne_plot()
                if (input$cb_cellranger) {
                    plot_dat <- get_cluster_dat_cellranger()
                } else {
                    plot_dat <- get_cluster_dat_seurat()
                }
                limits_x = range(plot_dat$tSNE_1)
                limits_y = range(plot_dat$tSNE_2)

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
                    scale_colour_gradient(low = 'lightgrey', high = 'blue') +
                    scale_x_continuous(limits = limits_x) +
                    scale_y_continuous(limits = limits_y) +
                    geom_point(size = 1) +
                    coord_fixed() +
                    theme_bw() +
                    theme(panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          legend.position = "none")

                gene_name = get_input_gene()
                gene_title = gene_name_title(gene_name,
                                                  species = dataset_info$species)
                data.use = data.frame(FetchData(
                    object = dataset_info$rdat,
                    vars.all = gene_name,
                ), check.names = FALSE)
                colnames(data.use) = gene_title
                ident.use = FetchData(dataset_info$rdat,
                                      vars.all = res_name)[,1]
                # plot_3 <- Seurat::VlnPlot(dataset_info$rdat,
                #                           group.by = res_name,
                #                           features.plot = get_input_gene(),
                #                           do.return = TRUE)
                plot_3 <- Seurat:::SingleVlnPlot(
                    feature = gene_title,
                    data = data.use,
                    cell.ident = ident.use,
                    gene.names = gene_title,
                    do.sort = FALSE,
                    y.max = NULL,
                    size.x.use = 16,
                    size.y.use = 16,
                    size.title.use = 20,
                    adjust.use = 1,
                    point.size.use = 1,
                    cols.use = NULL,
                    y.log = FALSE,
                    x.lab.rot = FALSE,
                    y.lab.rot = FALSE,
                    legend.position = 'right',
                    remove.legend = FALSE
                )

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
    }, height = plot_height_func(scale_factor_default))

    output$plot_gene_expr2 <- renderPlot({
        res_name = glue('res.{dataset_info$resolution}')
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

            plot_cor = cor(plot_dat[[get_input_gene1()]],
                           plot_dat[[get_input_gene2()]])
            plot_title = glue('{get_input_gene1()} & {get_input_gene2()} (cluster {paste(input$cluster_id, collapse = "/")}) [r = {round(plot_cor, digits = 2)}]')

            ggplot(plot_dat,
                aes_string(
                    x = get_input_gene1(),
                    y = get_input_gene2()
                    )
            ) +
            geom_point(size = 3, alpha = 0.7) +
            geom_rug(sides = 'bl') +
            ggtitle(plot_title) +
            coord_fixed() +
            theme_bw() +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  legend.position = "none")

        }
    }, height = plot_height_func(scale_factor_default - 0.1))

    find_marker_gene <- function() {
        withProgress(message = 'Find marker gene',
                     detail = 'collect cell id in cluster',
                     value = 0, {
            res_name = glue('res.{dataset_info$resolution}')
            cell_1 = rownames(dataset_info$rdat@meta.data)[dataset_info$rdat@meta.data[[res_name]] == input$sig_cluster_1]

            incProgress(amount = 0.2, message = 'run program')
            if (length(cell_1) <= 3) {
                runjs("document.getElementById('warning_info').innerHTML = '<font color=red>Warning: Cell group 1 has fewer than 3 cells</font>'")
                output_df = empty_sig_df
            } else {
                if (input$sig_cluster_2 == '(All other cells)') {
                    runjs("document.getElementById('warning_info').innerHTML = '<font color=blue>Everything OK!</font>'")
                    if (input$marker_pos == 'pos') {
                        output_df = FindMarkers(dataset_info$rdat,
                                                ident.1 = cell_1,
                                                ident.2 = NULL,
                                                test.use = 'roc',
                                                min.pct = 0.25,
                                                only.pos = TRUE
                        )
                    } else {
                        output_df = FindMarkers(dataset_info$rdat,
                                                ident.1 = cell_1,
                                                ident.2 = NULL,
                                                test.use = 'roc',
                                                min.pct = 0.25,
                                                only.pos = FALSE
                        )
                        output_df = subset(output_df, avg_diff < 0)
                    }
                } else {
                    cell_2 = rownames(dataset_info$rdat@meta.data)[dataset_info$rdat@meta.data[[res_name]] == input$sig_cluster_2]
                    if (length(cell_2) <= 3) {
                        runjs("document.getElementById('warning_info').innerHTML = '<font color=red>Warning: Cell group 2 has fewer than 3 cells</font>'")
                        output_df = empty_sig_df
                    } else {
                        runjs("document.getElementById('warning_info').innerHTML = '<font color=blue>Everything OK!</font>'")
                        if (input$marker_pos == 'pos') {

                            output_df = FindMarkers(dataset_info$rdat,
                                                    ident.1 = cell_1,
                                                    ident.2 = cell_2,
                                                    test.use = 'roc',
                                                    min.pct = 0.25,
                                                    only.pos = TRUE
                            )
                        } else {
                            runjs("document.getElementById('warning_info').innerHTML = '<font color=blue>Everything OK!</font>'")
                            output_df = FindMarkers(dataset_info$rdat,
                                                    ident.1 = cell_2,
                                                    ident.2 = cell_1,
                                                    test.use = 'roc',
                                                    min.pct = 0.25,
                                                    only.pos = TRUE
                            )
                        }
                    }
                }
            }

            incProgress(amount = 0.7, message = 'create output dataframe')

            if (!('gene' %in% colnames(output_df))) {
                output_df %<>%
                    rownames_to_column(var = 'gene') %>%
                    as_data_frame() %>%
                    dplyr::select(-p_val_adj) %>%
                    dplyr::filter(myAUC >= 0.7) %>%
                    arrange(desc(myAUC))
            }

            setProgress(value = 1, message = 'Finish!')
        })
        output_df
    }

    get_sig_cluster_input <- reactive({
        list(cluster1 = input$sig_cluster_1,
             cluster2 = input$sig_cluster_2)
    }) %>% debounce(2000)

    observeEvent({
        get_sig_cluster_input()
        input$marker_pos
        }, {
        if ((input$sig_cluster_1 != '') &&
            (input$sig_cluster_1 != input$sig_cluster_2)) {

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

    observe({
        query <- parseQueryString(session$clientData$url_search)
        if (!is.null(query[['dataset']])) {
            if (query[['dataset']] %in% resource_list$label) {
                updateSelectizeInput(session,
                                     inputId = 'dataset',
                                     selected = query[['dataset']])
            }
        }
    })
})

