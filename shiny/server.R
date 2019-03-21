# Sat Mar  9 23:31:25 2019 ------------------------------

library(shiny)
library(shinyjs)
library(tidyverse)
library(magrittr)
library(Seurat)
library(DT)
library(cowplot)
library(glue)
library(jsonlite)

resource.list <- read_json('data/resource_list.json')
source('config.txt')

# gene official symbols and their synonyms
# see: ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/
gene_human <- read_csv('gene/gene-human.csv', col_types = 'cc')
gene_mouse <- read_csv('gene/gene-mouse.csv', col_types = 'cc')

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

get_dr_name <- function(object, reduction.type = 'tsne') {
    dr.key <- GetDimReduction(
        object,
        reduction.type = reduction.type,
        slot = 'key'
    )
    return(paste0(dr.key, 1:2))
}

# Get data_frame for dimensional reduction
# here `object` is a Seurat object
get_dr_df <- function(object, reduction.type = 'tsne') {
    dr.df <- GetDimReduction(
        object,
        reduction.type = reduction.type,
        slot = 'cell.embeddings'
    ) %>%
        as.data.frame() %>%
        rownames_to_column('Barcode') %>%
        as_data_frame()

    dr.names <- get_dr_name(object, reduction.type)
    dr.df %>%
        rename_('dr_1' = dr.names[1], 'dr_2' = dr.names[2])
}


# here `object` is a Seurat object
get_cluster_df <- function(
    object,
    reduction.type = 'tsne',
    resolution.value = 0.8,
    object.parent = NULL,
    show.parent = TRUE,
    use.parent.dr = TRUE
) {
    if ((!is.null(object.parent)) && use.parent.dr) {
        dr.df <- get_dr_df(
            object.parent,
            reduction.type = reduction.type
        )
    } else {
        dr.df <- get_dr_df(
            object,
            reduction.type = reduction.type
        )
    }
    # dr.df <- get_dr_df(object, reduction.type = reduction.type)
    res.name = paste0('res.', resolution.value)
    cluster.df <- object@meta.data %>%
        rownames_to_column('Barcode') %>%
        as_data_frame() %>%
        dplyr::select(Barcode, one_of(res.name)) %>%
        dplyr::rename_(cluster = res.name) %>%
        dplyr::mutate(cluster = factor(cluster, levels = as.character(sort(as.integer(unique(cluster))))))

    cluster.dat <- left_join(
        dr.df, cluster.df, by = 'Barcode'
    )

    # print('pre-filter')
    # print(cluster.dat)
    if ((!is.null(object.parent)) && (!show.parent)) {
        cluster.dat <- cluster.dat %>%
            filter(!is.na(cluster))
    }

    # print('post-filter')
    # print(cluster.dat)
    cluster.dat
}

get_dataset_info <- function(object) {
    n_genes = nrow(object@data)
    n_cells = ncol(object@data)
    return(glue('{n_genes} genes across {n_cells} samples'))
}

sig.df.empty <- data_frame(
    gene = character(),
    myAUC = numeric(),
    avg_logFC = numeric(),
    pct.1 = numeric(),
    pct.2 = numeric()
)

shinyServer(function(input, output, session) {
    plot_height_func <- function(scale_factor) {
        return(function() {
            session$clientData$output_plot_gene_expr_width * scale_factor
        })
    }

    url_ctrl <- reactiveValues(
        # dataset = NULL,
        reso = NULL,
        dr = NULL,
        gene = NULL
    )

    dataset_info <- reactiveValues(
        resource = NULL,
        rdat = NULL,

        # subset data
        resource_subset = NULL,
        rdat_subset = NULL,

        # other tmp data
        img_name = '',
        info_text = '',
        deg_table = NULL
    )

    observeEvent(input$dataset, {
        if (input$dataset == 'none') {
            shinyjs::hide(id = 'dat_config')
            shinyjs::hide(id = 'dat_panel')

            dataset_info$resource = NULL
            dataset_info$rdat = NULL
            dataset_info$resource_subset = NULL
            dataset_info$rdat_subset = NULL
            dataset_info$info_text = ''
            dataset_info$img_name = ''
            dataset_info$deg_table = NULL

            updateSelectizeInput(session,
                                 inputId = 'sig_cluster_1',
                                 choices = NULL,
                                 selected = NULL)
            updateSelectizeInput(session,
                                 inputId = 'sig_cluster_2',
                                 choices = NULL,
                                 selected = NULL)
        } else {
            shinyjs::show(id = 'dat_config')
            shinyjs::show(id = 'dat_panel')
            shinyjs::hide(id = 'cb_allpt')
            shinyjs::hide(id = 'cluster_id_subset')
            shinyjs::hide(id = 'resolution_subset')

            message.main = 'Load seurat object from disk'
            withProgress(
                message = message.main,
                detail = 'Collect RDS file info',
                value = 0.1, {
                    resource.id = match(
                        input$dataset,
                        map_chr(resource.list, ~.$label)
                    )
                    resource = resource.list[[resource.id]]
                    if (is.null(resource$resolution)) {
                        resource$resolution = res_default
                    } else {
                        cluster.res = round(
                            resource$resolution,
                            digits = 1
                        )
                        if (cluster.res > 0 && cluster.res < 2 &&
                            cluster.res != 0.8 &&
                            is.null(url_ctrl$reso)) {
                            updateSliderInput(
                                session = session,
                                inputId = 'resolution',
                                value = cluster.res)
                        }
                    }

                    if (!is.null(url_ctrl$reso)) {
                        updateSliderInput(
                            session = session,
                            inputId = 'resolution',
                            value = url_ctrl$reso
                        )
                        url_ctrl$reso = NULL
                    }


                    if (is.null(resource$dims.use)) {
                        resource$dims.use = 1:dims_default
                    }
                    if (is.null(resource$species)) {
                        resource$species = species_default
                    }

                    if (!is.null(url_ctrl$dr)) {
                        updateSelectizeInput(
                            session = session,
                            inputId = 'dr_method',
                            selected = url_ctrl$dr
                        )
                        url_ctrl$dr = NULL
                    }

                    if (!is.null(url_ctrl$gene)) {
                        updateTextInput(
                            session = session,
                            inputId = 'tx_gene',
                            value = url_ctrl$gene
                        )
                        url_ctrl$gene = NULL
                    }

                    incProgress(
                        amount = 0.1,
                        message = 'Load seurat object',
                        detail = 'Read RDS file'
                    )

                    rdat = read_rds(
                        file.path(
                            'data',
                            input$dataset,
                            glue('{input$dataset}.rds')
                        )
                    )
                    incProgress(
                        amount = 0.6,
                        message = message.main,
                        detail = 'Get resolution list'
                    )

                    incProgress(
                        amount = 0.1,
                        message = message.main,
                        detail = 'Organize pre-defined subsets'
                    )
                    if (!is.null(resource$subset)) {
                        updateSelectizeInput(
                            session = session,
                            inputId = 'cb_subset',
                            choices = c(
                                '(None)' = 'none',
                                'Custom' = 'custom',
                                set_names(
                                    map_chr(resource$subset, ~.$label),
                                    map_chr(resource$subset, ~.$description)
                                )
                            )
                        )
                    }

                    setProgress(
                        value = 1,
                        message = message.main,
                        detail = 'All done!'
                    )
                })

            dataset_info$resource = resource
            dataset_info$rdat = rdat
            dataset_info$resource_subset = NULL
            dataset_info$rdat_subset = NULL
            dataset_info$img_name = ''
            dataset_info$info_text = get_dataset_info(rdat)
            dataset_info$deg_table = sig.df.empty
            runjs("document.getElementById('warning_info').innerHTML = ''")
        }
    })

    # ---------- 1st panel
    update_resolution_subset <- function(res) {
        res_name = glue('res.{res}')
        if (!is.null(dataset_info$rdat_subset) &&
            !(res_name %in% colnames(dataset_info$rdat_subset@meta.data))) {
            print(glue('update cluster resolution to {res}...'))
            message.main = 'Calculate clusters with new resolution'
            withProgress(
                message = message.main,
                detail = 'Run `FindClusters()` function',
                value = 0.1, {
                    dataset_info$rdat_subset = FindClusters(
                        object = dataset_info$rdat_subset,
                        dims.use = 1:dataset_info$resource_subset$dims.use,
                        print.output = FALSE,
                        resolution = res,
                        force.recalc = TRUE
                    )
                    setProgress(
                        value = 1,
                        message = message.main,
                        detail = 'All done!'
                    )
                }
            )
        }
    }

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
                                           dataset_info$resource$species,
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

    get_cluster_dat <- reactive({
        if (input$cb_subset != 'none' && (!is.null(dataset_info$rdat_subset))) {
            res_str = paste0('res.', input$resolution_subset)
            if (!(res_str %in%
                  colnames(dataset_info$rdat_subset@meta.data))) {
                dataset_info$rdat_subset = FindClusters(
                    object = dataset_info$rdat_subset,
                    dims.use = 1:dataset_info$resource$dims.use,
                    resolution = input$resolution_subset,
                    print.output = FALSE,
                    force.recalc = TRUE
                )
            }
            cluster.dat <- get_cluster_df(
                object = dataset_info$rdat_subset,
                reduction.type = input$dr_method,
                resolution.value = input$resolution_subset,
                object.parent = dataset_info$rdat,
                show.parent = input$cb_allpt
            )
        } else {
            cluster.dat <- get_cluster_df(
                object = dataset_info$rdat,
                reduction.type = input$dr_method,
                resolution.value = input$resolution,
            )

            # multiple selection
            updateSelectizeInput(
                session = session,
                inputId = 'cluster_id_subset',
                choices = levels(cluster.dat$cluster)
            )
        }

        updateSelectizeInput(
            session = session,
            inputId = 'cluster_id',
            choices = levels(cluster.dat$cluster)
        )

        # single selection
        updateSelectizeInput(
            session = session,
            inputId = 'sig_cluster_1',
            choices = c('', levels(cluster.dat$cluster)),
            selected = ''
        )
        # single selection
        updateSelectizeInput(
            session = session,
            inputId = 'sig_cluster_2',
            choices = c('(All other cells)', levels(cluster.dat$cluster)),
            selected = '(All other cells)'
        )

        cluster.dat
    })

    get_dim_plot <- reactive({

        plot_dat <- get_cluster_dat()
        # print(plot_dat)
        dr.names <- get_dr_name(
            dataset_info$rdat,
            reduction.type = input$dr_method
        )

        cls <- plot_dat$cluster
        levels(cls) <- paste0(
            levels(cls), ' (n=', table(cls)[levels(cls)], ')'
        )
        plot_dat$cluster_str = cls
        cluster.stats <- plot_dat %>%
            group_by(cluster) %>%
            summarise(
                dr_1 = median(dr_1),
                dr_2 = median(dr_2)
            )
        ggplot(data = plot_dat,
               mapping = aes(x = dr_1, y = dr_2, color = cluster_str)) +
            geom_point(size = 1) +
            scale_color_discrete(na.value = 'lightgrey') +
            geom_text(
                cluster.stats,
                mapping = aes(x = dr_1, y = dr_2, label = cluster),
                color = 'black',
                fontface = "bold",
                size = 6
            ) +
            coord_fixed() +
            labs(x = dr.names[1], y = dr.names[2]) +
            theme_bw() +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  legend.position = "none")
    })

    get_feature_plot <- reactive({
        req(get_input_gene())
        dr.names <- get_dr_name(
            dataset_info$rdat,
            reduction.type = input$dr_method
        )

        # limits_x = range(plot_dat$tSNE_1)
        # limits_y = range(plot_dat$tSNE_2)
        plot_dat <- left_join(
            get_cluster_dat(),
            enframe(dataset_info$rdat@data[get_input_gene(),],
                    name = 'Barcode', value = 'expr'),
            by = 'Barcode'
        )

        ggplot(plot_dat, aes(x = dr_1, y = dr_2, color = expr)) +
            scale_colour_gradient(low = 'lightgrey', high = 'blue') +
            geom_point(size = 1) +
            coord_fixed() +
            labs(x = dr.names[1], y = dr.names[2]) +
            theme_bw() +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  legend.position = "none")
    })

    get_vln_plot <- reactive({
        req(get_input_gene())

        gene_title = gene_name_title(
            get_input_gene(),
            species = dataset_info$resource$species
        )
        plot.dat <- left_join(
            get_cluster_dat(),
            enframe(dataset_info$rdat@data[get_input_gene(),],
                    name = 'Barcode', value = 'expr'),
            by = 'Barcode'
        )

        ggplot(plot.dat, aes(x = cluster, y = expr)) +
            geom_violin(aes(color = cluster, fill = cluster)) +
            stat_summary(
                geom = 'crossbar',
                fun.y = mean,
                fun.ymin = mean,
                fun.ymax = mean,
                width = 0.5
            ) +
            labs(x = 'Cluster ID', y = '', title = gene_title) +
            theme_bw() +
            theme(
                panel.grid = element_blank(),
                plot.title = element_text(hjust = 0.5, face = 'bold'),
                axis.title.x = element_blank(),
                axis.text.x = element_text(size = 16, face = 'bold')
            )
    })

    get_gene_expr_plot <- reactive({
        req(dataset_info$rdat)
        gene.name = get_input_gene()
        if (gene.name %in% rownames(dataset_info$rdat@data)) {
            plot_1 <- get_dim_plot()
            plot_2 <- get_feature_plot()
            plot_3 <- get_vln_plot()
            dataset_info$img_name = paste(
                dataset_info$resource$label,
                paste0('res', input$resolution),
                input$dr_method,
                gene.name,
                sep = '_'
            )

            plot_grid(
                plot_grid(plot_1, plot_2, align = 'h'),
                plot_3, ncol = 1,
                rel_heights = c(3, 2)
            )
        } else {
            dataset_info$img_name = paste(
                dataset_info$resource$label,
                paste0('res', input$resolution),
                input$dr_method,
                sep = '_'
            )
            get_dim_plot() +
                theme(
                    legend.position = 'right',
                    legend.title = element_blank(),
                    legend.text = element_text(size = 14),
                    axis.title = element_text(size = 16),
                    axis.text = element_text(size = 14)
                ) +
                guides(color = guide_legend(override.aes = list(size = 8)))
        }
    })

    output$plot_gene_expr <- renderPlot({
        get_gene_expr_plot()
        # }, width = 800, height = 600)
    }, height = plot_height_func(scale_factor_default))

    output$d_img <- downloadHandler(
        filename = function() {
            paste0(dataset_info$img_name, '.pdf')
        },
        content = function(file) {
            ggsave(filename = file, plot = get_gene_expr_plot(),
                   width = 9, height = 9)
        }
    )

    output$dat_info_text <- renderText({
        dataset_info$info_text
    })

    get_cluster_subset <- reactive({
        input$cluster_id_subset
    }) %>% debounce(2000)

    observe(if (is.null(input$cluster_id_subset)) {
        shinyjs::enable(id = 'resolution')

        dataset_info$rdat_subset = NULL
        dataset_info$resource_subset = NULL
        if (!is.null(dataset_info$rdat)) {
            dataset_info$info_text = get_dataset_info(dataset_info$rdat)
        }
    })

    observeEvent(get_cluster_subset(), {
        shinyjs::disable(id = 'resolution')

        message.main = 'Create new dataset using selected clusters'
        withProgress(
            message = message.main,
            detail = 'Make a copy of the original data',
            value = 0.1, {
                rdat_subset <- SetAllIdent(
                    dataset_info$rdat,
                    id = paste0('res.', input$resolution)
                )
                incProgress(
                    amount = 0.2,
                    message = message.main,
                    detail = 'Remove unused cells'
                )
                rdat_subset = SubsetData(
                    rdat_subset,
                    ident.use = get_cluster_subset(),
                    subset.raw = TRUE
                )
                incProgress(
                    amount = 0.3,
                    message = message.main,
                    detail = 'Clean invalid metadata'
                )
                rdat_subset@meta.data = rdat_subset@meta.data[,!str_detect(colnames(rdat_subset@meta.data), '^res\\.')]
                dataset_info$info_text = get_dataset_info(rdat_subset)
                incProgress(
                    amount = 0.1,
                    message = message.main,
                    detail = 'Find clusters on subset'
                )
                rdat_subset = FindClusters(
                    rdat_subset,
                    dims.use = 1:dataset_info$resource$dims.use,
                    resolution = res_default,
                    force.recalc = TRUE,
                    print.output = FALSE
                )
                setProgress(
                    value = 1,
                    message = message.main,
                    detail = 'All done!'
                )
            }
        )

        dataset_info$rdat_subset = rdat_subset
        resource_subset = list(
            species = dataset_info$resource$species,
            dims.use = dataset_info$resource$dims.use,
            resolution = dataset_info$resource$resolution
        )
        dataset_info$resource_subset = resource_subset
    })

    observeEvent(input$cb_subset, {
        if (input$cb_subset == 'none') {
            shinyjs::hide('cb_allpt')
            shinyjs::hide('cluster_id_subset')
            shinyjs::hide('resolution_subset')
            shinyjs::enable(id = 'resolution')
            dataset_info$rdat_subset = NULL
            dataset_info$resource_subset = NULL
            if (!is.null(dataset_info$rdat)) {
                dataset_info$info_text = get_dataset_info(dataset_info$rdat)
            }
            updateSelectizeInput(
                session = session,
                inputId = 'cluster_id_subset',
                selected = NULL
            )
            # updateCheckboxInput(
            #     session = session,
            #     inputId = 'cb_allpt',
            #     value = FALSE
            # )
            if (input$resolution_subset != res_default) {
                updateSliderInput(
                    session = session,
                    inputId = 'resolution_subset',
                    value = res_default
                )
            }
        } else if (input$cb_subset == 'custom') {
            shinyjs::show(id = 'cb_allpt')
            shinyjs::show(id = 'cluster_id_subset')
            shinyjs::show(id = 'resolution_subset')
        } else {
            # use pre-defined subsets
            shinyjs::show(id = 'cb_allpt')
            shinyjs::hide(id = 'cluster_id_subset')
            shinyjs::show(id = 'resolution_subset')

            message.main = 'Load seurat object from disk'
            withProgress(
                message = message.main,
                detail = 'Load RDS file',
                value = 0.1, {
                    rdat_subset = read_rds(
                        file.path(
                            'data',
                            input$dataset,
                            glue('{input$cb_subset}.rds')
                        )
                    )

                    incProgress(
                        amount = 0.6,
                        message = message.main,
                        detail = 'Collect dataset info'
                    )
                    dataset_info$rdat_subset = rdat_subset
                    dataset_info$info_text = get_dataset_info(rdat_subset)
                    resource_subset.id = match(
                        input$cb_subset,
                        map_chr(dataset_info$resource$subset, ~.$label)
                    )
                    resource_subset = dataset_info$resource$subset[[resource_subset.id]]
                    if (is.null(resource_subset$dims.use)) {
                        resource_subset$dims.use = 1:dataset_info$resource$dims.use
                    }
                    if (is.null(resource_subset$resolution)) {
                        resource_subset$resolution = res_default
                    }
                    dataset_info$resource_subset = resource_subset
                    setProgress(
                        value = 1,
                        message = message.main,
                        detail = 'All done!'
                    )
                }
            )
        }
    })

    # ---------- 2nd panel
    find_marker_gene <- function() {
        message.main = 'Find signature genes'
        withProgress(
            message = message.main,
            detail = 'Set cell ident using selected resolution',
            value = 0.1, {
            if (input$cb_subset != 'none' && (!is.null(dataset_info$rdat_subset))) {
                rdat <- SetAllIdent(
                    dataset_info$rdat_subset,
                    id = paste0('res.', input$resolution_subset)
                )
            } else {
                rdat <- SetAllIdent(
                    dataset_info$rdat,
                    id = paste0('res.', input$resolution)
                )
            }
            cell_1 = input$sig_cluster_1
            cell_2 = input$sig_cluster_2

            incProgress(
                amount = 0.2,
                message = message.main,
                detail = 'Finds markers for identity classes'
            )
            if (sum(rdat@ident %in% cell_1) <= 3) {
                runjs("document.getElementById('warning_info').innerHTML = '<font color=red>Warning: Cell group 1 has fewer than 3 cells</font>'")
                output_df = sig.df.empty
            } else {
                if (cell_2 == '(All other cells)') {
                    runjs("document.getElementById('warning_info').innerHTML = '<font color=blue>Everything OK!</font>'")
                    output_df = FindMarkers(
                        rdat,
                        ident.1 = cell_1,
                        test.use = 'roc',
                        min.pct = 0.5,
                        min.diff.pct = -Inf,
                        logfc.threshold = log(1.5),
                        only.pos = (input$marker_pos == 'pos'),
                        print.bar = FALSE
                    )
                    if (input$marker_pos != 'pos') {
                        output_df = subset(output_df, avg_diff < 0)
                    }
                } else {
                    if (sum(rdat@ident %in% cell_2) <= 3) {
                        runjs("document.getElementById('warning_info').innerHTML = '<font color=red>Warning: Cell group 2 has fewer than 3 cells</font>'")
                        output_df = empty_sig_df
                    } else {
                        runjs("document.getElementById('warning_info').innerHTML = '<font color=blue>Everything OK!</font>'")
                        output_df = FindMarkers(
                            rdat,
                            ident.1 = cell_1,
                            ident.2 = cell_2,
                            test.use = 'roc',
                            min.pct = 0.5,
                            min.diff.pct = -Inf,
                            logfc.threshold = log(1.5),
                            only.pos = (input$marker_pos == 'pos'),
                            print.bar = FALSE
                        )
                        if (input$marker_pos != 'pos') {
                            output_df = subset(output_df, avg_diff < 0)
                        }
                    }
                }
            }

            incProgress(
                amount = 0.6,
                message = message.main,
                detail = 'Create data frame for output'
            )

            output_df <- output_df %>%
                rownames_to_column(var = 'gene') %>%
                as_data_frame() %>%
                dplyr::select(
                    gene, myAUC, avg_logFC, pct.1, pct.2
                ) %>%
                dplyr::mutate(myAUC = pmax(myAUC, 1 - myAUC)) %>%
                dplyr::filter(myAUC > 0.7) %>%
                arrange(desc(myAUC))

            setProgress(
                value = 1,
                message = message.main,
                detail = 'All done!'
            )
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
            dataset_info$deg_table = find_marker_gene()
        }
    })

    output$table_sig_gene <- DT::renderDataTable({
        dataset_info$deg_table
    }, selection = 'single')

    observe(
        if (!is.null(input$table_sig_gene_row_last_clicked)) {
            gene_name = dataset_info$deg_table$gene[
                input$table_sig_gene_row_last_clicked]

            if (!is.na(gene_name)) {
                updateTextInput(session, "tx_gene", value = gene_name)
            }
        }
    )

    # ---------- 3rd panel
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
                                           dataset_info$resource$species,
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
                                           dataset_info$resource$species,
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

    output$plot_gene_expr2 <- renderPlot({
        req(dataset_info$rdat)
        if (input$cb_subset != 'none' && (!is.null(dataset_info$rdat_subset))) {
            rdat <- SetAllIdent(
                dataset_info$rdat_subset,
                id = paste0('res.', input$resolution_subset)
            )
        } else {
            rdat <- SetAllIdent(
                dataset_info$rdat,
                id = paste0('res.', input$resolution)
            )
        }

        if (get_input_gene1() != '' &&
            (get_input_gene1() %in% rownames(rdat@data)) &&
            get_input_gene2() != '' &&
            (get_input_gene2() %in% rownames(rdat@data)) &&
            get_input_gene1() != get_input_gene2()
        ) {
            valid_cells = rdat@cell.names[rdat@ident %in% input$cluster_id]
            plot_dat <- inner_join(
                enframe(
                    rdat@data[get_input_gene1(),],
                    name = 'Barcode', value = get_input_gene1()
                ),
                enframe(
                    rdat@data[get_input_gene2(),],
                    name = 'Barcode', value = get_input_gene2()
                ),
                by = 'Barcode'
            ) %>%
                dplyr::filter(Barcode %in% valid_cells)

            plot_dat <- left_join(
                plot_dat,
                enframe(rdat@ident, name = 'Barcode', value = 'cluster')
            )

            plot_cor = cor(plot_dat[[get_input_gene1()]],
                           plot_dat[[get_input_gene2()]])
            plot_title = glue('{get_input_gene1()} & {get_input_gene2()} (cluster {paste(input$cluster_id, collapse = "/")}) [r = {round(plot_cor, digits = 2)}]')

            ggplot(plot_dat,
                   aes_string(
                       x = get_input_gene1(),
                       y = get_input_gene2(),
                       color = 'cluster'
                   )
            ) +
                geom_point(size = 3, alpha = 0.7) +
                geom_rug(sides = 'bl') +
                ggtitle(plot_title) +
                coord_fixed() +
                theme_bw() +
                theme(
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    legend.title = element_blank(),
                    legend.position = "right",
                    legend.text = element_text(size = 14),
                    axis.title = element_text(size = 16),
                    axis.text = element_text(size = 14)
                ) +
                guides(color = guide_legend(override.aes = list(size = 8)))

        }
    }, height = plot_height_func(scale_factor_default - 0.1))

    # ---------- 4th panel
    output$plot_data_quality <- renderPlot({
        req(dataset_info$rdat)
        plot.dat <- left_join(
            get_cluster_dat(),
            dataset_info$rdat@meta.data[,c('nUMI', 'nGene')] %>%
                rownames_to_column('Barcode') %>%
                as_data_frame(),
            by = 'Barcode'
        )
        plot.umi <- ggplot(
            plot.dat,
            aes(x = cluster, y = nUMI)
        ) +
            geom_violin(aes(color = cluster, fill = cluster)) +
            stat_summary(
                geom = 'crossbar',
                fun.y = mean,
                fun.ymin = mean,
                fun.ymax = mean,
                width = 0.5
            ) +
            labs(x = '', y = 'nUMI', title = 'nUMI') +
            theme_bw() +
            theme(
                panel.grid = element_blank(),
                legend.position = 'none',
                axis.text = element_text(size = 14),
                axis.title = element_text(size = 16),
                plot.title = element_text(
                    size = 18, face = 'bold', hjust = 0.5
                )
            )
        plot.gene <- ggplot(
            plot.dat,
            aes(x = cluster, y = nGene)
        ) +
            geom_violin(aes(color = cluster, fill = cluster)) +
            stat_summary(
                geom = 'crossbar',
                fun.y = mean,
                fun.ymin = mean,
                fun.ymax = mean,
                width = 0.5
            ) +
            labs(x = 'Cluster ID', y = 'nGene', title = 'nGene') +
            theme_bw() +
            theme(
                panel.grid = element_blank(),
                legend.position = 'none',
                axis.text = element_text(size = 14),
                axis.title = element_text(size = 16),
                plot.title = element_text(
                    size = 18, face = 'bold', hjust = 0.5
                )
            )
        plot_grid(plot.umi, plot.gene, ncol = 1)
    })

    # ---------- parsing URL
    # /?dataset=cb&dr=tsne&reso=0.3&gene=AXL
    observe({
        query <- parseQueryString(session$clientData$url_search)
        if (!is.null(query[['dataset']])) {
            if (query[['dataset']] %in% map_chr(resource.list, ~.$label)) {
                updateSelectizeInput(
                    session,
                    inputId = 'dataset',
                    selected = query[['dataset']])

                if (!is.null(query[['dr']]) &&
                    (query[['dr']] %in% c('pca', 'tsne', 'umap'))) {
                    url_ctrl$dr = query[['dr']]
                }

                reso <- as.numeric(query[['reso']])
                if (length(reso) > 0 && !is.na(reso) &&
                    reso >= 0 && reso <= 2) {
                    url_ctrl$reso = round(reso, digits = 1)
                }

                if (is.character(query[['gene']])) {
                    url_ctrl$gene = query[['gene']]
                }
            }
        }
    })
})
