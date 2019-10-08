########################################################
# DETERMING OPTIMAL NETWORK MODULES FROM EACH COARSE CELL TYPE
# VERSION 2.
#
########################################################

library(WGCNA)
library(tibble)
library(ape)
library(factoextra)
library(glmnet)
library(FNN)
library(dplyr)
library(cluster)
library(MASS)

###############################################################
# LOADING THE DATA
# I need some sort of way to read in what the cell types are...or just type them
coarseCellTypes = c('CD4.T.cells', 'CD8.T.cells', 'NK.cells', 'B.cells', 'monocytic.lineage', 'neutrophils', 'endothelial.cells', 'fibroblasts')

# this is annoying...
for(ct in coarseCellTypes){
  assign(paste(ct), read.csv(paste0('../mixed_dfs/',ct,'.csv'), row.names = 1))
}
# so i'm doing it manually and saving the data frames
fibroblasts <- read.csv('../mixed_dfs/fibroblasts.csv', row.names = 1)
NK.cells <- t(NK.cells)

cell_df_list = list('bcell'=B.cells, "cd4"=CD4.T.cells, "cd8"=CD8.T.cells, 'endo'=endothelial.cells,
                    "fibro"=fibroblasts, "mono"=monocytic.lineage, "neutro"=neutrophils, "nkcells"=NK.cells)
save(cell_df_list, file='coarse_cell_data.RData')

load("coarse_cell_data.RData")

###############################################################

#####
# Creating and saving an object with gene module identity and eigengenes...
#
# WORKING WITH A SINGLE COARSE CELL TYPE
#   [1] "b"      "CD4"    "CD8"    "nk"     "neutro" "mono"   "fibro"  "endo"  
#
# First trial with neutrophils. They seem to be the most different in tSNE
#####
powers = c(c(1:10), seq(from = 12, to=20, by=2))

s = abs(bicor(cell_df_list[['nkcells']]))
sft = pickSoftThreshold(cell_df_list[['nkcells']], powerVector = powers, verbose = 5)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab='Soft Threshold (power)',ylab='Scale Free Topology Model Fit,signed R^2',
     type='n', main = paste('Scale independence'))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=1,col='red'); abline(h=0.90,col='red')

# correlations look strange...
beta = 12
a = s^beta
rm(s)
# w = 1-a #this is an alternative to TOM
TOM = TOMsimilarity(a)
rm(a)
dissTOM = 1-TOM
rm(TOM)
geneTree = hclust(as.dist(dissTOM), method = 'average')
modules = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 4, pamRespectsDendro = FALSE,
                        minClusterSize = 30)
module.colours = labels2colors(modules)
plotDendroAndColors(geneTree, module.colours, 'Module colours', dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, main='')
MEs = moduleEigengenes(cell_df_list[['nkcells']], colors = module.colours, excludeGrey = FALSE)$eigengenes
MEDiss = 1-cor(MEs)
METree = hclust(as.dist(MEDiss), method = "average")
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

MEDissThres = 0.02
abline(h=MEDissThres, col = "red")
merge = mergeCloseModules(cell_df_list[['nkcells']], module.colours, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
plotDendroAndColors(geneTree, cbind(module.colours, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

nkcells_colors = mergedColors
nkcells_MEs = mergedMEs
rm(dissTOM)

cell_mod_colors = list('bcell' = bcell_colors, 'cd4' = cd4_colors, 'cd8' = cd8_colors, 'endo' = endo_colors,
                       'neutro' = neutro_colors, 'fibro'=fibro_colors, 'mono'=mono_colors, 'nkcells'=nkcells_colors)
save(cell_mod_colors, file = 'coarse_cell_mod_colors.RData')
# THIS PROCESS WAS REPEATED FOR ALL COARSE CELL TYPES

#####
# CALCULATING MODULE ENRICHMENT...
#   1. make a list of all modules
#   2. How to choose the best module for a cell type?
#     -sum module gene probabilities assuming hypergeometric distribution?
#     -use module size?
#     -lowest covariance between module genes and other cell types?
#     -fold change expression of each gene in the module to other cell types?
#     -lowest heirarchical clustering distance?
#     -best silhouette score
#   3. save genes of module
#
#####

#####
# testing 2. silhouette score 
#####

# function to get joined gene symbol list, sort based on reference, impute expression of any missing gene symbols
prep_df <- function(ct_df, ref_df){
  # browser()
  ref_genes = colnames(ref_df)
  test_genes = colnames(ct_df)
  common_genes = intersect(ref_genes, test_genes)
  # drop irrelevant genes
  ct_df <- ct_df[,common_genes]
  # add new genes
  # though this may not even apply to some comparisons...
  new_genes <- ref_genes[!(ref_genes %in% common_genes)]
  print(length(new_genes))
  # impute
  # OPTION 1: impute values from expected distance from sample median
  #   this is really slow...
  # print('making impute difference')
  # ref_gene_diff <- sapply(new_genes, function(gs) mean(apply(ref_df, 1, function(gsm) median(gsm)-gsm[gs])))
  # print('making new data frame')
  # new_df <- sapply(new_genes, function(gs) apply(ct_df, 1, function(gsm) median(gsm) - ref_gene_diff[gs]))
  
  # OPTION 2: just use the median expression
  # new_df <- sapply(new_genes, function(x) apply(ct_df, 1, median))
  if(length(new_genes) > 0){
    median_values <- apply(ct_df, 1, median)
    new_df <- matrix(rep(median_values, length(new_genes)), ncol = length(new_genes))
    colnames(new_df) <- new_genes
    # sort
    ct_df <- cbind(ct_df, new_df)
  }
  ct_df <- ct_df[,colnames(ref_df)]
  return(ct_df)
}

# in a test of 'prep_df' almost 10% of the genes had to be imputed...this could lead to false negatives
# [1] "bcell"   "cd4"     "cd8"     "endo"    "fibro"   "mono"    "neutro"  "nkcells"
test_df_list <- cell_df_list[!(names(cell_df_list) %in% 'cd4')]
target_df <- cell_df_list[['cd4']]

# load the diluted data frames
file_pattern <- sapply(seq(0,75,25), function(x) paste0('cd4tcells_coarse_',x,'.csv'))
diluted_dfs <- lapply(file_pattern, function(f) t(read.csv(paste0('../mixed_dfs/',f), row.names = 1)))
names(diluted_dfs) <- seq(0,75,25)

# make all the dfs have the same genes as target df
adjusted_df_list <- list()
for(i in 1:length(test_df_list)){
  ct_df <- prep_df(test_df_list[[i]], target_df)
  adjusted_df_list[[names(test_df_list)[i]]] <- ct_df
}

# data frame to store silhouette results
mod_res <- data.frame(matrix(NA, ncol = 6, nrow = n_distinct(cell_mod_colors[['cd4']])))
names(mod_res) <- c('module', 'mod_size', 'ct_score', 'ct_avg_score', 'dilution_score', 'd_avg_score')
mod_res$module <- unique(cell_mod_colors[['cd4']])
counter <- 1

# module scores against other cell types:
for(mod_color in unique(cell_mod_colors[['cd4']])){
  print(paste('progress:', counter/n_distinct(cell_mod_colors[['cd4']])))
  target_mod_genes <- cell_mod_colors[['cd4']] == mod_color
  mod_res[mod_res$module == mod_color, 'mod_size'] <- sum(target_mod_genes)
  target_mod_df <- t(target_df[, target_mod_genes])
  test_mod_dfs <- lapply(adjusted_df_list, function(df) t(df[, target_mod_genes]))
  all_mod_dfs <- append(test_mod_dfs, list('cd4' = target_mod_df), 0)
  clusterCodes <- vector()
  for(i in 1:length(all_mod_dfs)){
    df_length <- rep(i, ncol(all_mod_dfs[[i]]))
    clusterCodes <- c(clusterCodes, df_length)
  }
  all_ct_df <- do.call(cbind, all_mod_dfs)
  # browser()
  mod_silhouette <- silhouette(clusterCodes, dist(t(all_ct_df)))
  s <- summary(mod_silhouette)
  # s_value <- mean(s$clus.avg.widths)/var(s$clus.avg.widths)
  mod_res[mod_res$module == mod_color, 'ct_score'] <- s$clus.avg.widths[1]
  mod_res[mod_res$module == mod_color, 'ct_avg_score'] <- s$avg.width
  counter <- counter + 1
}

# module scores against diluted cell type

adjusted_df_list <- list()
for(i in 1:length(diluted_dfs)){
  ct_df <- prep_df(diluted_dfs[[i]], target_df)
  adjusted_df_list[[names(diluted_dfs)[i]]] <- ct_df
}

counter <- 1
for(mod_color in unique(cell_mod_colors[['cd4']])){
  print(paste('progress:', counter/n_distinct(cell_mod_colors[['cd4']])))
  target_mod_genes <- cell_mod_colors[['cd4']] == mod_color
  # mod_res[mod_res$module == mod_color, 'mod_size'] <- sum(target_mod_genes)
  target_mod_df <- t(target_df[, target_mod_genes])
  test_mod_dfs <- lapply(adjusted_df_list, function(df) t(df[, target_mod_genes]))
  all_mod_dfs <- append(test_mod_dfs, list('cd4' = target_mod_df), 0)
  clusterCodes <- vector()
  for(i in 1:length(all_mod_dfs)){
    df_length <- rep(i, ncol(all_mod_dfs[[i]]))
    clusterCodes <- c(clusterCodes, df_length)
  }
  all_ct_df <- do.call(cbind, all_mod_dfs)
  # browser()
  mod_silhouette <- silhouette(clusterCodes, dist(t(all_ct_df)))
  s <- summary(mod_silhouette)
  # s_value <- mean(s$clus.avg.widths)/var(s$clus.avg.widths)
  mod_res[mod_res$module == mod_color, 'dilution_score'] <- s$clus.avg.widths[1]
  mod_res[mod_res$module == mod_color, 'd_avg_score'] <- s$avg.width
  counter <- counter + 1
}

# 3...
b_mod_geneNames = colnames(coarse_df)[b_colors=='pink']

#####
# 1. CREATING ME GLMNET MODELS
# 2. CREATING MEAN TOPOLOGY DISTANCES
#
#####

# 1...
b_pink_df = b_cell_df[,b_colors=='pink']
b_ME_model = cv.glmnet(as.matrix(b_pink_df), matrix(b_MEs$MEpink, nrow = nrow(b_pink_df), ncol = 1), 
                          family = 'gaussian', alpha = 0.5)

# 2...
# change: ct_ (11), and ME color (7)

b_ME_dist = matrix(0, nrow = ncol(b_pink_df), ncol = ncol(b_pink_df))
for(i in 1:ncol(b_pink_df)){
  gene_ME_dist = matrix(0, nrow = ncol(b_pink_df), ncol = nrow(b_pink_df))
  for(j in 1:nrow(b_pink_df)){
    gene_scaler = abs(b_pink_df[j,i]-b_MEs$MEpink[j])
    ME_dist = sapply(unlist(b_pink_df[j,]), function(x) abs(x - b_MEs$MEpink[j])/gene_scaler)
    gene_ME_dist[,j] = ME_dist
  }
  b_avg_dist = apply(gene_ME_dist, 1, mean)
  b_ME_dist[i,] <- b_avg_dist
}


par(mfrow=c(3,2))
for(s in sample(1:72,6)){
  hist(b_ME_dist[,s], breaks=20)
}

apply(b_ME_dist, 2, mean)
apply(b_ME_dist, 2, var)

x = b_ME_dist[,60]
h <- hist(x, breaks = 20, freq = F)
x.fit <- seq(min(x), max(x), length=length(x))
# plot(h)
# fitdistr(x, 'Poisson')
# curve(dpois(x, 0.655), from = 0, to = 10) #strange looking
# fitdistr(x, 'Beta')
# curve(dbeta(x, ))
fitdistr(x, 'Gamma')
x.gamma <- dgamma(x.fit, shape=0.77, rate=0.805)
lines(x.fit, x.gamma/max(x.gamma)*max(h$density), col='green')
# curve(dgamma(x, shape=1.275, rate=1.944), from = 0, to = 10) #looks good
fitdistr(x, 'Weibull')
x.weibull <- dweibull(x.fit, shape=0.797, scale=0.814)
lines(x.fit, x.weibull/max(x.weibull)*max(h$density), col='red')
# curve(dweibull(x, shape=0.978, scale=0.647), from = 0, to = 10) #looks okay...

# b_ME_dist_norm = t(apply(b_ME_dist, 1, function(x) {x/(var(x))}))
# b_ME_dist_norm = apply(b_ME_dist_norm, 2, function(x) {x/(mean(x))})
# 
# apply(b_ME_dist_norm, 1, mean)
# apply(b_ME_dist_norm, 1, var)
# 
# hist(b_ME_dist_norm[,45], breaks=20)

# b_samp_avg_dist = apply(b_ME_dist, 2, mean)
# hist(b_samp_avg_dist, breaks = 20)
b_gene_avg_dist = apply(b_ME_dist, 1, mean)
hist(b_gene_avg_dist, breaks = 20)

hist(b_avg_dist, breaks=20)




# for predicting a single sample
# predict(mono_ME_model, s=mono_ME_model$lambda.1se, newx = t(as.matrix(unlist(mono_steelblue_df[3,]))))

save(b_avg_dist, cd4_avg_dist, cd8_avg_dist, endo_avg_dist, fibro_avg_dist, mono_avg_dist, neutro_avg_dist, nk_avg_dist,
     b_mod_geneNames, cd4_mod_geneNames, cd8_mod_geneNames, endo_mod_geneNames, fibro_mod_geneNames, mono_mod_geneNames,
     neutro_mod_geneNames, nk_mod_geneNames, b_ME_model, cd4_ME_model, cd8_ME_model, endo_ME_model, fibro_ME_model,
     mono_ME_model, neutro_ME_model, nk_ME_model, file='WGCNA_module_output.RData')

