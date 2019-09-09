# FUNCTIONS USED IN 'decon_prediction.R'

format_data <- function(test_df, scale){
  # browser()
  print("data frame dimensions:")
  print(dim(test_df))
  #scale test_df
  if(scale == 'Linear'){
    test_df = log2(test_df)
  }else if(scale == 'Log10'){
    test_df = test_df*log2(10)
  }
  #impute NAs
  print("# of NA's:")
  num.na <- sum(is.na(test_df))
  print(num.na)
  if(num.na > 0){
    na.rows <- unname(which(apply(test_df, 1, function(x) sum(is.na(x))==ncol(test_df))))
    test_df <- test_df[-na.rows,]
    na.impute <- unname(which(apply(test_df, 1, function(x) sum(is.na(x))>0)))
    for(r in na.impute){
      test_df[r,is.na(test_df[r,])] <- median(unlist(test_df[r,]), na.rm = T)
    }
  }
  #format
  test_df = t(test_df)
  mix_cols = colnames(test_df)[colnames(test_df) %in% colnames(coarse_df)]
  test_df = test_df[,mix_cols]
  missing_genes = colnames(coarse_df)[(!(colnames(coarse_df) %in% colnames(test_df)))]
  median_expression = apply(test_df,1,median)
  missing_genes_df = sapply(missing_genes, function(g) median_expression, simplify = 'array')
  test_df_all = cbind(test_df, missing_genes_df)
  test_df_all = test_df_all[,colnames(coarse_df)]
  return(test_df_all)
}

make_prediction <- function(cell_model, test_df){
  # browser()
  test_df_mod = test_df[,colnames(cell_model$df)]
  test_dist = array(matrix(0, nrow = ncol(test_df_mod), ncol = ncol(test_df_mod)),
                      dim = c(ncol(test_df_mod),ncol(test_df_mod),nrow(test_df_mod)))
  for(i in 1:nrow(test_df_mod)){
    gene_ME_dist = matrix(0, nrow = ncol(test_df_mod), ncol = ncol(test_df_mod))
    for(j in 1:ncol(test_df_mod)){
      ME_dist = sapply(unlist(test_df_mod[i,]), function(x) (x-test_df_mod[i,j])^2)
      gene_ME_dist[j,] = ME_dist
    }
    test_dist[,,i] <- gene_ME_dist
  }
  
  test_rotation = t(sapply(1:dim(test_dist)[3], function(i) svd(test_dist[,,i])$u[,1]))
  
  test_pred <- predict(cell_model$model, s=cell_model$model$lambda.1se, newx=test_rotation)
  # print(quantile(test_pred))
  # print('done')
  return(test_pred)
}





