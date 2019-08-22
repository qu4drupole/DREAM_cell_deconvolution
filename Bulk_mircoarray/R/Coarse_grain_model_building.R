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
#   [1] "b"      "CD4"    "CD8"    "nk"     "neutro" "mono"   "fibro"  "endo"  
#
#####

cd4_cell_df = coarse_df[cellTypes=='CD4.T.cells',]
cd8_cell_df = coarse_df[cellTypes=='CD8.T.cells',]
mono_cell_df = coarse_df[cellTypes=='monocytic.lineage',]
endo_cell_df = coarse_df[cellTypes=='endothelial.cells',]
fibro_cell_df = coarse_df[cellTypes=='fibroblasts',]
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
  nk_mod_test = t(nk_cell_df[,neutro_colors==mod_color])
  cd4_mod_test = t(cd4_cell_df[,neutro_colors==mod_color])
  cd8_mod_test = t(cd8_cell_df[,neutro_colors==mod_color])
  mono_mod_test = t(mono_cell_df[,neutro_colors==mod_color])
  endo_mod_test = t(endo_cell_df[,neutro_colors==mod_color])
  fibro_mod_test = t(fibro_cell_df[,neutro_colors==mod_color])
  clusterCodes = c(rep(1,ncol(neutro_mod_test)), rep(2,ncol(b_mod_test)), rep(3,ncol(nk_mod_test)),
                   rep(4,ncol(cd4_mod_test)), rep(5,ncol(cd8_mod_test)), rep(6,ncol(mono_mod_test)),
                   rep(7,ncol(endo_mod_test)), rep(8,ncol(fibro_mod_test)))
  coarse_mod_df = cbind(neutro_mod_test, b_mod_test, nk_mod_test, cd4_mod_test, cd8_mod_test, 
                        mono_mod_test, endo_mod_test, fibro_mod_test)
  mod_silhouette <- silhouette(clusterCodes, dist(t(coarse_mod_df)))
  s = summary(mod_silhouette)
  print(mean(s$clus.avg.widths)/var(s$clus.avg.widths))
}

# NOTES:
# CD8: I chose a slightly lower score in favor a larger module size (~50 vs. 150)
# fibro: I chose second highest score; max had 800 genes
# nk: a lot of the options kind of sucked...the final module is quite large (740)

# c/p -- cellType

neutro_mod_df = neutro_cell_df[,neutro_colors=='orange']

neutro_dist = matrix(0, nrow = ncol(neutro_mod_df), ncol = ncol(neutro_mod_df))
for(i in 1:ncol(neutro_mod_df)){
  gene_dist_matrix = matrix(0, nrow = ncol(neutro_mod_df), ncol = nrow(neutro_mod_df))
  for(j in 1:nrow(neutro_mod_df)){
    gene_dist = sapply(unlist(neutro_mod_df[j,]), function(x) x/neutro_mod_df[j,i])
    gene_dist_matrix[,j] = gene_dist
  }
  neutro_avg_dist = apply(gene_dist_matrix, 1, mean)
  neutro_dist[i,] <- neutro_avg_dist
}

save(cd4_mod_df, cd4_dist, cd8_mod_df, cd8_dist, mono_mod_df, mono_dist, endo_mod_df, endo_dist, 
     fibro_mod_df, fibro_dist, b_mod_df, b_dist, nk_mod_df, nk_dist, neutro_mod_df, neutro_dist, 
     file='coarse_cell_training_data.RData')


#####
# PREPARING THE IN SILICO CELL MIXTURES
#
#####

# c/p -- cellType, mix_int
# !! CAUTION: neutro file types have different naming convention !!
# !! CAUTION: all mix files have different naming conventions !!

#//// TO DO ////#
# Something is wrong with the 0% b cell file--remake
#//// TO DO ////#

nk_30mix = read.csv('../mixed_dfs/nkcells_coarse_all_30.csv')
nk_30mix <- as.data.frame(nk_30mix)
row.names(nk_30mix) = nk_30mix$gene_symbol_sql
nk_30mix = subset(nk_30mix, select = -c(gene_symbol_sql))
nk_30mix = t(nk_30mix)
mix_cols = colnames(nk_30mix)[colnames(nk_30mix) %in% colnames(coarse_df)]
nk_30mix = nk_30mix[,mix_cols]
missing_genes = colnames(coarse_df)[(!(colnames(coarse_df) %in% colnames(nk_30mix)))]
median_expression = apply(nk_30mix,1,median)
missing_genes_df = sapply(missing_genes, function(g) median_expression, simplify = 'array')
nk_30mix_all = cbind(nk_30mix, missing_genes_df)
nk_30mix_all = nk_30mix_all[,colnames(coarse_df)]
nk_30mix_mod = nk_30mix_all[,colnames(nk_mod_df)]

nk_30_dist = array(matrix(0, nrow = ncol(nk_30mix_mod), ncol = ncol(nk_30mix_mod)),
                         dim = c(ncol(nk_30mix_mod),ncol(nk_30mix_mod),nrow(nk_30mix_mod)))
for(i in 1:nrow(nk_30mix_mod)){
  gene_ME_dist = matrix(0, nrow = ncol(nk_30mix_mod), ncol = ncol(nk_30mix_mod))
  for(j in 1:ncol(nk_30mix_mod)){
    ME_dist = sapply(unlist(nk_30mix_mod[i,]), function(x) x/nk_30mix_mod[i,j])
    gene_ME_dist[j,] = ME_dist
  }
  nk_30_dist[,,i] <- gene_ME_dist
}

#\/\/\/\/ PROBLEM \/\/\/\/#
#
# "mono" has a lot of 'NA' values
#
#\/\/\/\/ PROBLEM \/\/\/\/#


# c/p -- cellType

fibroArrayList = list('100' = fibro_100_dist, '90' = fibro_100_dist, '80' = fibro_80_dist,
                       '70' = fibro_70_dist, '60' = fibro_60_dist, '50' = fibro_50_dist,
                       '40' = fibro_40_dist, '30' = fibro_30_dist, '20' = fibro_20_dist,
                       '10' = fibro_10_dist, '0' = fibro_0_dist)


initial = 0
for(idx in 1:length(fibroArrayList)){
  if(initial == 0){
    all_cov_matrix = t(sapply(1:100, function(i) sapply(1:ncol(fibro_mod_df), 
                                                        function(x) cov(fibro_dist[x,],fibroArrayList[[idx]][x,,i]))))
    all_cov_matrix = cbind(all_cov_matrix, rep(as.integer(names(fibroArrayList)[idx]),100))
    colnames(all_cov_matrix) = c(colnames(fibro_mod_df), "mix")
    initial = 1
  } else {
    cov_matrix = t(sapply(1:100, function(i) sapply(1:ncol(fibro_mod_df), 
                                                    function(x) cov(fibro_dist[x,],fibroArrayList[[idx]][x,,i]))))
    cov_matrix = cbind(cov_matrix, rep(as.integer(names(fibroArrayList)[idx]),100))
    colnames(cov_matrix) = c(colnames(fibro_mod_df), "mix")
    all_cov_matrix = rbind(all_cov_matrix, cov_matrix)
  }
}

#//// CAUTION ////#
#
# "endo" had a few 'NA's that I replaced with matrix mean
# "fibro" had a few 'NA's that I replaced with matrix mean
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

fibro_model <- fit10

# NOTES:
# fibro sucks...

#############
# --STOPPING POINT--
# DIDN'T FINISH "nk" OR 'neutro'
# PROBABLY NEED TO PICK A DIFFERENT NK MODULE--ONE MIX IS 0.5 Gb
#############

save(cd4_mod_df, cd4_dist, cd8_mod_df, cd8_dist, mono_mod_df, mono_dist, endo_mod_df, endo_dist, 
     fibro_mod_df, fibro_dist, b_mod_df, b_dist, nk_mod_df, nk_dist, neutro_mod_df, neutro_dist,
     cd4_model, cd8_model, mono_model, endo_model, fibro_model, b_model, nk_model, neutro_model,
     file='coarse_cell_data_models.RData')
