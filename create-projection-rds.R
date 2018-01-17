#!/usr/bin/env Rscript
# Thu Jan 18 04:06:58 2018 ------------------------------
library(tidyverse)

projection_file = commandArgs(TRUE)[1]
output_rds_file = commandArgs(TRUE)[2]
read_csv(projection_file) %>%
    mutate(Barcode = str_extract(Barcode, '^[^-]+')) %>%
    rename(tSNE_1 = `TSNE-1`, tSNE_2 = `TSNE-2`) %>%
    write_rds(output_rds_file)
