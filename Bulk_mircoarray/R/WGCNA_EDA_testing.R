########################################################
# DETERMING OPTIMAL NETWORK MODULES FROM EACH COARSE CELL TYPE
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
# Creating and saving an object with gene module identity and eigengenes...
#
# WORKING WITH A SINGLE COARSE CELL TYPE
#   [1] "b"      "CD4"    "CD8"    "nk"     "neutro" "mono"   "fibro"  "endo"  
#
# First trial with neutrophils. They seem to be the most different in tSNE
#####
powers = c(c(1:10), seq(from = 12, to=20, by=2))

neutro_cell_df = coarse_df[cellTypes=='neutrophils',]
s = abs(bicor(neutro_cell_df))
sft = pickSoftThreshold(neutro_cell_df, powerVector = powers, verbose = 5)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab='Soft Threshold (power)',ylab='Scale Free Topology Model Fit,signed R^2',
     type='n', main = paste('Scale independence'))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=1,col='red'); abline(h=0.90,col='red')

beta = 12
a = s^beta
# w = 1-a
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
neutro_MEs = mergedMEs

# THIS PROCESS WAS REPEATED FOR ALL COARSE CELL TYPES

#####
# CALCULATING MODULE ENRICHMENT...
#   1. make a list of all modules
#   2. How to choose the best module for a cell type?
#     Xsum module gene probabilities assuming hypergeometric distribution?
#     -use module size?
#     -lowest covariance between module genes and other cell types?
#     -fold change expression of each gene in the module to other cell types?
#     -lowest heirarchical clustering distance?
#     -best silhouette score
#     -k-medoids
#   3. save genes of module
#
#####

#####
# testing 2. silhouette score with neutro vs. b cells and NK cells
#####
cd4_cell_df = coarse_df[cellTypes=='CD4.T.cells',]
cd8_cell_df = coarse_df[cellTypes=='CD8.T.cells',]
mono_cell_df = coarse_df[cellTypes=='monocytic.lineage',]
endo_cell_df = coarse_df[cellTypes=='endothelial.cells',]
fibro_cell_df = coarse_df[cellTypes=='fibroblasts',]
b_cell_df = coarse_df[cellTypes=='B.cells',]
nk_cell_df = coarse_df[cellTypes=='NK.cells',]

# calculating the Frobenius norm which is essentially the euclidean norm...?
norm(neutro_mod_test, "F") #627.96
norm(b_mod_test, "F") #1129.926
norm(nk_mod_test, "F") #328.81

#testing all neutro modules
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
  print(mean(s$clus.avg.widths))
}

# I think the higher the score the more unique the silhouette 
# Based on the scores and module sizes (shooting for ~100), "salmon" is my top choice

# 'darkturquoise' is the best

#####
# 1. CREATING ME GLMNET MODELS
# 2. CREATING MEAN TOPOLOGY DISTANCES
#
#####

# 1...
neutro_darkturquoise_df = neutro_cell_df[,neutro_colors=='darkturquoise']
neutro_ME_model = cv.glmnet(as.matrix(neutro_darkturquoise_df), matrix(neutro_MEs$MEdarkturquoise, nrow = nrow(neutro_darkturquoise_df), ncol = 1), 
                          family = 'gaussian', alpha = 0.5)

# how does the ME work for the other cell types?
b_darkturquoise_df = b_cell_df[,neutro_colors=='darkturquoise']
b_darkturquoise_ME = predict(neutro_ME_model, s=neutro_ME_model$lambda.1se, newx = as.matrix(b_darkturquoise_df))

nk_darkturquoise_df = nk_cell_df[,neutro_colors=='darkturquoise']
nk_darkturquoise_ME = predict(neutro_ME_model, s=neutro_ME_model$lambda.1se, newx = as.matrix(nk_darkturquoise_df))

#pick some random genes
par(mfrow=c(3,3))
for(g in sample(1:sum(neutro_colors=='darkturquoise'),3)){
  plot(neutro_MEs$MEdarkturquoise~neutro_darkturquoise_df[,g], main="neutro")
  plot(b_darkturquoise_ME~b_darkturquoise_df[,g], main="b cells")
  plot(nk_darkturquoise_ME~nk_darkturquoise_df[,g], main="nk cells")
}

# testing the module on in silico mixtures

neutro_100mix = read.csv('../mixed_dfs/neutro_coarse_100.csv')
# a little restructring...
neutro_100mix <- as.data.frame(neutro_100mix)
row.names(neutro_100mix) = neutro_100mix$gene_symbol_sql
neutro_100mix = subset(neutro_100mix, select = -c(gene_symbol_sql))
neutro_100mix = t(neutro_100mix)
# a little more formatting
mix_cols = colnames(neutro_100mix)[colnames(neutro_100mix) %in% colnames(coarse_df)]
neutro_100mix = neutro_100mix[,mix_cols]
missing_genes = colnames(coarse_df)[(!(colnames(coarse_df) %in% colnames(neutro_100mix)))]
median_expression = apply(neutro_100mix,1,median)
missing_genes_df = sapply(missing_genes, function(g) median_expression, simplify = 'array')
neutro_100mix_all = cbind(neutro_100mix, missing_genes_df)
neutro_100mix_all = neutro_100mix_all[,colnames(coarse_df)]
neutro_100mix_mod = neutro_100mix_all[,neutro_colors=='darkturquoise']
neutro_100mix_ME = predict(neutro_ME_model, s=neutro_ME_model$lambda.1se, newx = as.matrix(neutro_100mix_mod))

par(mfrow=c(3,5))
for(g in sample(1:sum(neutro_colors=='darkturquoise'),3)){
  plot(neutro_MEs$MEdarkturquoise~neutro_darkturquoise_df[,g], main="training")
  plot(neutro_100mix_ME~neutro_100mix_mod[,g], main="100% neutro")
  plot(neutro_75mix_ME~neutro_75mix_mod[,g], main="75% neutro")
  plot(neutro_50mix_ME~neutro_50mix_mod[,g], main="50% neutro")
  plot(neutro_25mix_ME~neutro_25mix_mod[,g], main="25% neutro")
}

#some genes separate convolutions well...others don't

# 2...
# repeat loop for each mixture

neutro25_ME_dist = matrix(0, nrow = ncol(neutro_25mix_mod), ncol = ncol(neutro_25mix_mod))
for(i in 1:ncol(neutro_25mix_mod)){
  gene_ME_dist = matrix(0, nrow = ncol(neutro_25mix_mod), ncol = nrow(neutro_25mix_mod))
  for(j in 1:nrow(neutro_25mix_mod)){
    gene_scaler = abs(neutro_25mix_mod[j,i]-neutro_25mix_ME[j])
    ME_dist = sapply(unlist(neutro_25mix_mod[j,]), function(x) abs(x - neutro_25mix_ME[j])/gene_scaler)
    gene_ME_dist[,j] = ME_dist
  }
  neutro_avg_dist = apply(gene_ME_dist, 1, mean)
  neutro25_ME_dist[i,] <- neutro_avg_dist
}


par(mfrow=c(5,3))
par(mar=c(1,1,1,1))
neutroMixList = list('training' = neutro_ME_dist, '100' = neutro100_ME_dist, '75' = neutro75_ME_dist,
                     '50' = neutro50_ME_dist, '25' = neutro25_ME_dist)

randomGenes <- sample(1:sum(neutro_colors=='darkturquoise'),3)
for(data in 1:length(neutroMixList)){
  for(g in randomGenes){
    x = unname(neutroMixList[[data]])[,g]
    x.fit <- seq(min(x), max(x), length=length(x))
    h = hist(neutroMixList[[data]][,g], breaks=20, xlim = c(0.5,2), main = names(neutroMixList[data]))
    fit = fitdistr(neutroMixList[[data]][,g], 'Gamma')
    x.gamma <- dgamma(x.fit, shape=fit$estimate[1], rate=fit$estimate[2])
    lines(x.fit, x.gamma/max(x.gamma)*max(h$density), col='green', lwd=2)
  }
}

#####
# 1. Calculate a divergence/entropy score from fitted distributions
# 
#
#####

KL.divergence(neutro_ME_dist[,10],neutro100_ME_dist[,10],20)




# for predicting a single sample
# predict(mono_ME_model, s=mono_ME_model$lambda.1se, newx = t(as.matrix(unlist(mono_steelblue_df[3,]))))

save(b_avg_dist, cd4_avg_dist, cd8_avg_dist, endo_avg_dist, fibro_avg_dist, mono_avg_dist, neutro_avg_dist, nk_avg_dist,
     b_mod_geneNames, cd4_mod_geneNames, cd8_mod_geneNames, endo_mod_geneNames, fibro_mod_geneNames, mono_mod_geneNames,
     neutro_mod_geneNames, nk_mod_geneNames, b_ME_model, cd4_ME_model, cd8_ME_model, endo_ME_model, fibro_ME_model,
     mono_ME_model, neutro_ME_model, nk_ME_model, file='WGCNA_module_output.RData')

