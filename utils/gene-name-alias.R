library(tidyverse)
library(glue)

load('~/Project/gene-trend/shiny/gene/robj/human.RData')
load('~/Project/gene-trend/shiny/gene/robj/mouse.RData')

gene_info.h %>%
    dplyr::select(gene = Symbol, alias = Synonyms) %>%
    # dplyr::mutate(title = map2_chr(gene, alias, function(a, b) {
    #     if (b == '-') return(a)
    #     else {
    #         glue('{a} ({b})')
    #         return(glue('{a} ({b})'))
    #     }
    # })) %>%
    dplyr::filter(alias != '-') %>%
    dplyr::arrange(gene) %>%
    write_csv('shiny/gene/gene-human.csv')

gene_info.m %>%
    dplyr::select(gene = Symbol, alias = Synonyms) %>%
    # dplyr::mutate(title = map2_chr(gene, alias, function(a, b) {
    #     if (b == '-') return(a)
    #     else {
    #         glue('{a} ({b})')
    #         return(glue('{a} ({b})'))
    #     }
    # })) %>%
    dplyr::filter(alias != '-') %>%
    dplyr::arrange(gene) %>%
    write_csv('shiny/gene/gene-mouse.csv')
