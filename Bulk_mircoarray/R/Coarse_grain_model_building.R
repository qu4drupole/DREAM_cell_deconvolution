########################################################
# BUILDING THE REQUIRED DATA AND MODELS FOR DECONVOLUTION
# OF ALL COARSE GRAIN CELL TYPES
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

coarse_data = read.csv('RMA_coarse_081419.csv')
cellTypes = coarse_data$cellType
coarse_data = subset(coarse_data, select = -c(cellType))
coarse_df = as.data.frame(coarse_data)
rownames(coarse_df) = coarse_df$X
coarse_df = subset(coarse_df, select = -c(X))

#####
# Calculating the WGCNA modules for each cell type
#
# WORKING WITH A SINGLE COARSE CELL TYPE
#   [1] "b"      "CD4"    "CD8"    "neutro"     "neutro" "mono"   "b"  "endo"  
#
#####

cd4_cell_df = coarse_df[cellTypes=='CD4.T.cells',]
cd8_cell_df = coarse_df[cellTypes=='CD8.T.cells',]
mono_cell_df = coarse_df[cellTypes=='monocytic.lineage',]
endo_cell_df = coarse_df[cellTypes=='endothelial.cells',]
b_cell_df = coarse_df[cellTypes=='bblasts',]
b_cell_df = coarse_df[cellTypes=='B.cells',]
nk_cell_df = coarse_df[cellTypes=='NK.cells',]
neutro_cell_df = coarse_df[cellTypes=='neutrophils',]

# c/p -- cellType
powers = c(c(1:10), seq(from = 12, to=20, by=2))

s = abs(bicor(neutro_cell_df))
sft = pickSoftThreshold(neutro_cell_df, powerVector = powers, verbose = 5)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab='Soft Threshold (power)',ylab='Scale Free Topology Model Fit,signed R^2',
     type='n', main = paste('Scale independence'))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=1,col='red'); abline(h=0.90,col='red')

beta = 12
a = s^beta
TOM = TOMsimilarity(a)
dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = 'average')
modules = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 4, pamRespectsDendro = FALSE,
                        minClusterSize = 30)
module.colours = labels2colors(modules)
plotDendroAndColors(geneTree, module.colours, 'Module colours', dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, main='')
MEs = moduleEigengenes(neutro_cell_df, colors = module.colours, excludeGrey = FALSE)$eigengenes
MEDiss = 1-cor(MEs)
METree = hclust(as.dist(MEDiss), method = "average")
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.05
abline(h=MEDissThres, col = "red")
merge = mergeCloseModules(neutro_cell_df, module.colours, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
plotDendroAndColors(geneTree, cbind(module.colours, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

neutro_colors = mergedColors

#####
# CALCULATING MODULE SILHOUETTE SCORES 
#   
#####

# c/p -- cellType_colors
for(mod_color in unique(neutro_colors)){
  print(mod_color)
  print(sum(neutro_colors==mod_color))
  neutro_mod_test = t(neutro_cell_df[,neutro_colors==mod_color])
  b_mod_test = t(b_cell_df[,neutro_colors==mod_color])
  neutro_mod_test = t(neutro_cell_df[,neutro_colors==mod_color])
  cd4_mod_test = t(cd4_cell_df[,neutro_colors==mod_color])
  cd8_mod_test = t(cd8_cell_df[,neutro_colors==mod_color])
  mono_mod_test = t(mono_cell_df[,neutro_colors==mod_color])
  endo_mod_test = t(endo_cell_df[,neutro_colors==mod_color])
  b_mod_test = t(b_cell_df[,neutro_colors==mod_color])
  clusterCodes = c(rep(1,ncol(neutro_mod_test)), rep(2,ncol(b_mod_test)), rep(3,ncol(neutro_mod_test)),
                   rep(4,ncol(cd4_mod_test)), rep(5,ncol(cd8_mod_test)), rep(6,ncol(mono_mod_test)),
                   rep(7,ncol(endo_mod_test)), rep(8,ncol(b_mod_test)))
  coarse_mod_df = cbind(neutro_mod_test, b_mod_test, neutro_mod_test, cd4_mod_test, cd8_mod_test, 
                        mono_mod_test, endo_mod_test, b_mod_test)
  mod_silhouette <- silhouette(clusterCodes, dist(t(coarse_mod_df)))
  s = summary(mod_silhouette)
  print(mean(s$clus.avg.widths)/var(s$clus.avg.widths))
}

# NOTES:
# CD8: I chose a slightly lower score in favor a larger module size (~50 vs. 150)
# b: I chose second highest score; max had 800 genes
# neutro: a lot of the options kind of sucked...the final module is quite large (740)

# c/p -- cellType

neutro_mod_df = neutro_cell_df[,neutro_colors=='orange']

###
# If there are memory issues with the module, trim:
nk_mod_df = nk_mod_df[,sample(1:740,300)]
###

nk_dist = matrix(0, nrow = ncol(nk_mod_df), ncol = ncol(nk_mod_df))
for(i in 1:ncol(nk_mod_df)){
  gene_dist_matrix = matrix(0, nrow = ncol(nk_mod_df), ncol = nrow(nk_mod_df))
  for(j in 1:nrow(nk_mod_df)){
    gene_dist = sapply(unlist(nk_mod_df[j,]), function(x) (x-nk_mod_df[j,i])^2)
    gene_dist_matrix[,j] = gene_dist
  }
  nk_avg_dist = apply(gene_dist_matrix, 1, mean)
  nk_dist[i,] <- nk_avg_dist
}

save(cd4_mod_df, cd4_dist, cd8_mod_df, cd8_dist, mono_mod_df, mono_dist, endo_mod_df, endo_dist, 
     b_mod_df, b_dist, b_mod_df, b_dist, neutro_mod_df, neutro_dist, neutro_mod_df, neutro_dist,
     cellTypes, coarse_df, nk_mod_df, nk_dist, fibro_mod_df, file='coarse_cell_training_data.RData')



#####
# PREPARING THE IN SILICO CELL MIXTURES
#
#####

# c/p -- cellType, mix_int
# !! CAUTION: neutro file types have different naming convention !!
# !! CAUTION: all mix files have different naming conventions !!

#//// TO DO ////#
# Something is wrong with the 0% b cell file--remake
# Same thing seems to be true with neutro...maybe the mixer just fucks up with all 0's
#//// TO DO ////#

#replace 'cell.type' name
neutroArrayList <- list()
for(mixRatio in seq(0,100,10)){
  print(paste('working on',as.character(mixRatio)))
  ct_Nmix = read.csv(gsub('N',as.character(mixRatio), '../mixed_dfs/neutro_coarse_all_N.csv'))
  ct_Nmix <- as.data.frame(ct_Nmix)
  row.names(ct_Nmix) = ct_Nmix$gene_symbol_sql
  ct_Nmix = subset(ct_Nmix, select = -c(gene_symbol_sql))
  ct_Nmix = t(ct_Nmix)
  mix_cols = colnames(ct_Nmix)[colnames(ct_Nmix) %in% colnames(coarse_df)]
  ct_Nmix = ct_Nmix[,mix_cols]
  missing_genes = colnames(coarse_df)[(!(colnames(coarse_df) %in% colnames(ct_Nmix)))]
  median_expression = apply(ct_Nmix,1,median)
  missing_genes_df = sapply(missing_genes, function(g) median_expression, simplify = 'array')
  ct_Nmix_all = cbind(ct_Nmix, missing_genes_df)
  ct_Nmix_all = ct_Nmix_all[,colnames(coarse_df)]
  ct_Nmix_mod = ct_Nmix_all[,colnames(neutro_mod_df)]
  
  ct_N_dist = array(matrix(0, nrow = ncol(ct_Nmix_mod), ncol = ncol(ct_Nmix_mod)),
                           dim = c(ncol(ct_Nmix_mod),ncol(ct_Nmix_mod),nrow(ct_Nmix_mod)))
  for(i in 1:nrow(ct_Nmix_mod)){
    gene_ME_dist = matrix(0, nrow = ncol(ct_Nmix_mod), ncol = ncol(ct_Nmix_mod))
    for(j in 1:ncol(ct_Nmix_mod)){
      ME_dist = sapply(unlist(ct_Nmix_mod[i,]), function(x) (x-ct_Nmix_mod[i,j])^3) #took out square
      gene_ME_dist[j,] = ME_dist
    }
    ct_N_dist[,,i] <- gene_ME_dist
  }
  print(sum(is.na(ct_N_dist))/length(ct_N_dist))
  print(sum(is.infinite(ct_N_dist))/length(ct_N_dist))
  ct_N_list <- list(ct_N_dist)
  names(ct_N_list) <- as.character(mixRatio)
  neutroArrayList <- c(neutroArrayList, ct_N_list)
}



#\/\/\/\/ PROBLEM \/\/\/\/#
#
# Fibro had some NA's in 100%
#   -there must be a 0 in the expression data, which creates Inf of NaN in the ratio distance metric
fibroArrayList[['100']][is.na(fibroArrayList[['100']])] <- mean(fibroArrayList[['100']], na.rm=T)
fibroArrayList[['100']][is.infinite(fibroArrayList[['100']])] <- 0
#
#\/\/\/\/ PROBLEM \/\/\/\/#


# c/p -- cellType
# 
# fibroArrayList = list('100' = fibro_100_dist, '90' = fibro_100_dist, '80' = fibro_80_dist,
#                        '60' = fibro_60_dist, '50' = fibro_50_dist, '70' = fibro_70_dist,
#                        '40' = fibro_40_dist, '30' = fibro_30_dist, '20' = fibro_20_dist,
#                        '10' = fibro_10_dist, '0' = fibro_0_dist)


# 
# par(mfrow=c(3,2))
# par(mar=c(1,1,1,1))
# for(i in seq(1,11, by=2)){
#   # mean.dist <- vector()
#   # for(x in 1:ncol(neutro_mod_df)){
#   #   diff.gene <- (neutro_dist[x,] - neutroArrayList[[i]][x,,4])**2
#   #   log.gene <- log(diff.gene)
#   #   log.gene[log.gene==-Inf] <- 0
#   #   mean.dist <- append(mean.dist, exp(sum(log.gene))/ncol(neutro_dist))
#   # }
#   mean.dist <- sapply(1:ncol(neutro_mod_df), function(x) mean((neutro_dist[,x] - neutroArrayList[[i]][,x,7])**2)*100)
#   hist(mean.dist, breaks=10, xlim = c(0,20), main = names(neutroArrayList[[i]]))
# }
# 
# 
# exp(sum(log((neutro_dist[x,] - neutro_100_dist[x,,1])**2))/ncol(neutro_dist))
# exp(sum(log((neutro_dist[x,] - neutro_20_dist[x,,1])**2))/ncol(neutro_dist))
# exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
##################### 


initial = 0
for(idx in 1:length(neutroArrayList)){
  if(initial == 0){
    all_cov_matrix = t(sapply(1:100, function(i) svd(neutroArrayList[[idx]][,,i])$u[,1]))
    all_cov_matrix = cbind(all_cov_matrix, rep(as.integer(names(neutroArrayList)[idx]),100))
    colnames(all_cov_matrix) = c(colnames(neutro_mod_df), "mix")
    initial = 1
  } else {
    cov_matrix = t(sapply(1:100, function(i) svd(neutroArrayList[[idx]][,,i])$u[,1]))
    cov_matrix = cbind(cov_matrix, rep(as.integer(names(neutroArrayList)[idx]),100))
    colnames(cov_matrix) = c(colnames(neutro_mod_df), "mix")
    all_cov_matrix = rbind(all_cov_matrix, cov_matrix)
  }
}

#//// CAUTION ////#
#
# "cd8" skipping 70% mix because it had too many missing values
sum(is.na(all_cov_matrix))/length(all_cov_matrix)
all_cov_matrix[is.na(all_cov_matrix)] <- mean(all_cov_matrix, na.rm = T)
#
#//// CAUTION ////#

test_idx <- sample(1:1000,100)
test_cov_matrix <- all_cov_matrix[test_idx,]
train_cov_matrix <- all_cov_matrix[-test_idx,]



#####
# CREATING THE GLMNET MODELS FROM IN SILICO MIXTURES
#
#####

for (i in 0:10) {
  assign(paste("fit", i, sep=""), cv.glmnet(subset(train_cov_matrix, select = -c(mix)), train_cov_matrix[,'mix'], 
                                            type.measure="mse", alpha=i/10,family="gaussian"))
}

yhat0 <- predict(fit0, s=fit0$lambda.1se, newx=subset(test_cov_matrix, select = -c(mix)))
yhat1 <- predict(fit1, s=fit1$lambda.1se, newx=subset(test_cov_matrix, select = -c(mix)))
yhat2 <- predict(fit2, s=fit2$lambda.1se, newx=subset(test_cov_matrix, select = -c(mix)))
yhat3 <- predict(fit3, s=fit3$lambda.1se, newx=subset(test_cov_matrix, select = -c(mix)))
yhat4 <- predict(fit4, s=fit4$lambda.1se, newx=subset(test_cov_matrix, select = -c(mix)))
yhat5 <- predict(fit5, s=fit5$lambda.1se, newx=subset(test_cov_matrix, select = -c(mix)))
yhat6 <- predict(fit6, s=fit6$lambda.1se, newx=subset(test_cov_matrix, select = -c(mix)))
yhat7 <- predict(fit7, s=fit7$lambda.1se, newx=subset(test_cov_matrix, select = -c(mix)))
yhat8 <- predict(fit8, s=fit8$lambda.1se, newx=subset(test_cov_matrix, select = -c(mix)))
yhat9 <- predict(fit9, s=fit9$lambda.1se, newx=subset(test_cov_matrix, select = -c(mix)))
yhat10 <- predict(fit10, s=fit10$lambda.1se, newx=subset(test_cov_matrix, select = -c(mix)))

mean((test_cov_matrix[,"mix"] - yhat0)^2)
mean((test_cov_matrix[,"mix"] - yhat1)^2)
mean((test_cov_matrix[,"mix"] - yhat2)^2)
mean((test_cov_matrix[,"mix"] - yhat3)^2)
mean((test_cov_matrix[,"mix"] - yhat4)^2)
mean((test_cov_matrix[,"mix"] - yhat5)^2)
mean((test_cov_matrix[,"mix"] - yhat6)^2)
mean((test_cov_matrix[,"mix"] - yhat7)^2)
mean((test_cov_matrix[,"mix"] - yhat8)^2)
mean((test_cov_matrix[,"mix"] - yhat9)^2)
mean((test_cov_matrix[,"mix"] - yhat10)^2)

# c/p -- cellType, fitN
test_sample = sample(1:100,5)
test_cov_matrix[,"mix"][test_sample]
yhat4[test_sample]
plot(test_cov_matrix[,"mix"], yhat5)

fibro_model <- fit0


#############
# FINISHED: 
#   -neutro
#   -cd4 
#   -cd8 
#     -- skipping 70% mix because it had too many missing values
#   -mono
#     -- just a bad fit for the admix
#   -endo
#   -b
#     -- just a bad fit for the admix
#   -fibro
#     -- squared difference seems to work better than squared ratio
#     -- how about cubed difference? Including negative numbers really messes it up...
#   -nk
#
# --STOPPING POINT--
# 
#############

save(cd4_mod_df, cd4_model, neutro_mod_df, neutro_model, cd8_mod_df, cd8_model, mono_mod_df, mono_model,
     endo_mod_df, endo_model, fibro_mod_df, fibro_model, b_mod_df, b_model, nk_mod_df, nk_model, coarse_df,
     file='coarse_cell_data_models.RData')

cd4 = list('df' = cd4_mod_df, 'model' = cd4_model)
cd8 = list('df' = cd8_mod_df, 'model' = cd8_model)
neutro = list('df' = neutro_mod_df, 'model' = neutro_model)
mono = list('df' = mono_mod_df, 'model' = mono_model)
endo = list('df' = endo_mod_df, 'model' = endo_model)
fibro = list('df' = fibro_mod_df, 'model' = fibro_model)
b = list('df' = b_mod_df, 'model' = b_model)
nk = list('df' = nk_mod_df, 'model' = nk_model)

cell_models = list('CD4.T.cells' = cd4, 'CD8.T.cells' = cd8, 'NK.cells' = nk,
                   'B.cells' = b, 'monocytic.lineage' = mono, 'neutrophils' = neutro,
                   'endothelial.cells' = endo, 'fibroblasts' = fibro)

save(cell_models, coarse_df, file='coarse_cell_models_list.RData')









  
