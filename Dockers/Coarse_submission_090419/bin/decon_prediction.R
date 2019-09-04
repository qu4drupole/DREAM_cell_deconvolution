####################################################
# PRIMARY SCRIPT TO RUN DECONVOLUTION PREDICTIONS
#
# REQUIRES: 
#   'coarse_cell_data_models.RData' -- cell models and dfs
#   'functions.R' -- processing and formatting functions
#
####################################################

library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(tibble)
library(reshape2)
library(glmnet)

source("functions.R")
load("coarse_cell_models_list.RData")

# read meta data
input_df <- readr::read_csv("input/input.csv")
dataset_names <- input_df$dataset.name
expression_files  <- input_df$hugo.expr.file
expression_paths <- paste0("input/", expression_files)
data_scales <- input_df$scale

# need to nest module df and cell model in list
# e.g. cell_model = list('df'=cd4_mod_df, 'model' = cd4_model)
expression_path <- expression_paths[1]
dataset_name <- dataset_names[1]
data_scale <- data_scales[1]


do_decon <- function(expression_path, dataset_name, data_scale){
  #format the data
  ex <- expression_path %>%
    readr::read_csv() %>%
    as.data.frame() %>%
    tibble::column_to_rownames("Gene") %>%
    format_data(scale = data_scale)
  
  #calculate distance, rotation, and make predictions
  pred_matrix <- sapply(cell_models, function(x) make_prediction(x, ex))
  pred_df <- as.data.frame(pred_matrix)
  pred_df <- mutate(pred_df, 'sample.id' = rownames(ex), 'dataset.name' = dataset_name)
  pred_df <- melt(pred_df, id.vars = c('dataset.name', 'sample.id'), variable.name = 'cell.type', 
                  value.name = 'measured') %>%
    arrange(sample.id)
  return(pred_df)
}

test_res <- pmap(list(expression_paths, dataset_names, data_scales), do_decon)

combined_result_df <- dplyr::bind_rows(test_res)

dir.create("output")
readr::write_csv(combined_result_df, "output/predictions.csv")

