# Docker test script
 

library(readr)
library(tibble)
# library(WGCNA)

# print(list.files())
# print(getwd())
input_df <- read_csv("input/input.csv") #deleted readr::

## Extract the names of each dataset
dataset_names <- input_df$dataset.name

## Extract the names of the expression files that use 
## Hugo symbols as gene identifiers
expression_files  <- input_df$hugo.expr.file

## Form the paths of the expression files
expression_paths <- paste0("input/", expression_files)

cell_types = c("B.cells", "CD4.T.cells", "CD8.T.cells", "NK.cells", "neutrophils", "monocytic.lineage", "fibroblasts", "endothelial.cells")

ds_counter <- 0
for(f in expression_files){
  eds <- read_csv(paste("input/",f,sep=""))
  eds <- as.data.frame(eds)
  out_nrow <- length(cell_types) * sum(grepl('Sample', colnames(eds)))
  dataset_name = unlist(strsplit(f,".", fixed = T))[1]
  sample_names = rep(grep("Sample", colnames(eds), value=T),1,each=length(cell_types))
  cell.type.fill = rep(cell_types, length(grep("Sample", colnames(eds))))
  tmp_df <- data.frame("dataset.name" = dataset_name, 
                       "sample.id" = sample_names,
                       "cell.type" = cell.type.fill, 
                       "prediction" = numeric(out_nrow))
  if(ds_counter == 0){
    for(s in sample_names){
      for(ct in cell_types){
        if(ct == 'neutrophils'){
          tmp_df[(tmp_df$sample.id == s & tmp_df$cell.type == ct),'prediction'] <- 69
        } else {
          tmp_df[(tmp_df$sample.id == s & tmp_df$cell.type == ct),'prediction'] <- mean(sample(eds[,which(colnames(eds) == s)],10))
        }
      }
      out_df <- tmp_df
      ds_counter <- ds_counter + 1
    }
  } else {
    for(s in sample_names){
      for(ct in cell_types){
        if(ct == 'neutrophils'){
          tmp_df[(tmp_df$sample.id == s & tmp_df$cell.type == ct),'prediction'] <- 69
        } else {
          tmp_df[(tmp_df$sample.id == s & tmp_df$cell.type == ct),'prediction'] <- mean(sample(eds[,which(colnames(eds) == s)],10))
        }
      }
    }
    out_df = rbind(out_df, tmp_df)
  }
}

dir.create("output")
write_csv(out_df, "output/predictions.csv") # deleted "readr::"



