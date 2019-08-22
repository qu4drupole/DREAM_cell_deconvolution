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
# library(LaplacesDemon)

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
  print(mean(s$clus.avg.widths)/var(s$clus.avg.widths))
}

# I think the higher the score the more unique the silhouette 
# Based on the scores and module sizes (shooting for ~100), "salmon" is my top choice

# 'darkturquoise' and 'orange' are both good
neutro_darkturquoise_df = neutro_cell_df[,neutro_colors=='darkturquoise']
neutro_orange_df = neutro_cell_df[,neutro_colors=='orange']

# making the 'darkturquoise' gene topology
neutro_ME_dist = matrix(0, nrow = ncol(neutro_orange_df), ncol = ncol(neutro_orange_df))
for(i in 1:ncol(neutro_orange_df)){
  # gene_ME_dist = matrix(0, nrow = ncol(neutro_orange_df), ncol = nrow(neutro_orange_df))
  gene_dist_matrix = matrix(0, nrow = ncol(neutro_orange_df), ncol = nrow(neutro_orange_df))
  for(j in 1:nrow(neutro_orange_df)){
    # gene_scaler = neutro_orange_df[j,i]-neutro_MEs$MEorange[j]
    # ME_dist = sapply(unlist(neutro_orange_df[j,]), function(x) ((x - neutro_MEs$MEorange[j])/gene_scaler)**2)
    gene_dist = sapply(unlist(neutro_orange_df[j,]), function(x) x/neutro_orange_df[j,i])
    gene_dist_matrix[,j] = gene_dist
  }
  neutro_avg_dist = apply(gene_dist_matrix, 1, mean)
  neutro_ME_dist[i,] <- neutro_avg_dist
}

# Flipping the above loop to get gene-specific distances per sample
# neutro_ME_dist = matrix(0, nrow = nrow(neutro_darkturquoise_df), ncol = ncol(neutro_darkturquoise_df))
# for(i in 1:nrow(neutro_darkturquoise_df)){
#   # browser()
#   gene_ME_dist = matrix(0, nrow = ncol(neutro_darkturquoise_df), ncol = ncol(neutro_darkturquoise_df))
#   for(j in 1:ncol(neutro_darkturquoise_df)){
#     gene_scaler = abs(neutro_darkturquoise_df[i,j]-neutro_MEs$MEdarkturquoise[i])
#     ME_dist = sapply(unlist(neutro_darkturquoise_df[i,]), function(x) abs(x - neutro_MEs$MEdarkturquoise[i])/gene_scaler)
#     gene_ME_dist[,j] = ME_dist
#   }
#   neutro_avg_dist = apply(gene_ME_dist, 1, mean)
#   neutro_ME_dist[i,] <- neutro_avg_dist
# }

#####
# 1. CREATING ME GLMNET MODELS
# 2. CREATING MEAN TOPOLOGY DISTANCES
#
#####

# 1...

neutro_ME_model = cv.glmnet(as.matrix(neutro_orange_df), matrix(neutro_MEs$MEorange, nrow = nrow(neutro_orange_df), ncol = 1), 
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

neutro_90mix = read.csv('../mixed_dfs/neutro_coarse_all_90.csv')
# a little restructring...
neutro_90mix <- as.data.frame(neutro_90mix)
row.names(neutro_90mix) = neutro_90mix$gene_symbol_sql
neutro_90mix = subset(neutro_90mix, select = -c(gene_symbol_sql))
neutro_90mix = t(neutro_90mix)
# a little more formatting
mix_cols = colnames(neutro_90mix)[colnames(neutro_90mix) %in% colnames(coarse_df)]
neutro_90mix = neutro_90mix[,mix_cols]
missing_genes = colnames(coarse_df)[(!(colnames(coarse_df) %in% colnames(neutro_90mix)))]
median_expression = apply(neutro_90mix,1,median)
missing_genes_df = sapply(missing_genes, function(g) median_expression, simplify = 'array')
neutro_90mix_all = cbind(neutro_90mix, missing_genes_df)
neutro_90mix_all = neutro_90mix_all[,colnames(coarse_df)]
neutro_90mix_mod = neutro_90mix_all[,neutro_colors=='darkturquoise']
neutro_90mix_ME = predict(neutro_ME_model, s=neutro_ME_model$lambda.1se, newx = as.matrix(neutro_90mix_mod))

# redoing some stuff for the orange module

# neutro_90mix_mod = neutro_90mix_all[,neutro_colors=='orange']
# neutro_70mix_mod = neutro_70mix_all[,neutro_colors=='orange']
# neutro_50mix_mod = neutro_50mix_all[,neutro_colors=='orange']
# neutro_30mix_mod = neutro_30mix_all[,neutro_colors=='orange']
# neutro_10mix_mod = neutro_10mix_all[,neutro_colors=='orange']
# neutro_100mix_ME = predict(neutro_ME_model, s=neutro_ME_model$lambda.1se, newx = as.matrix(neutro_100mix_mod))
# neutro_80mix_ME = predict(neutro_ME_model, s=neutro_ME_model$lambda.1se, newx = as.matrix(neutro_80mix_mod))
# neutro_60mix_ME = predict(neutro_ME_model, s=neutro_ME_model$lambda.1se, newx = as.matrix(neutro_60mix_mod))
# neutro_40mix_ME = predict(neutro_ME_model, s=neutro_ME_model$lambda.1se, newx = as.matrix(neutro_40mix_mod))
# neutro_20mix_ME = predict(neutro_ME_model, s=neutro_ME_model$lambda.1se, newx = as.matrix(neutro_20mix_mod))

par(mfrow=c(3,6))
par(mar=c(1,1,1,1))
for(g in sample(1:sum(neutro_colors=='orange'),3)){
  plot(neutro_MEs$MEorange~neutro_orange_df[,g], main="training")
  plot(neutro_100mix_ME~neutro_100mix_mod[,g], main="100% neutro")
  plot(neutro_80mix_ME~neutro_80mix_mod[,g], main="80% neutro")
  plot(neutro_60mix_ME~neutro_60mix_mod[,g], main="60% neutro")
  plot(neutro_40mix_ME~neutro_40mix_mod[,g], main="40% neutro")
  plot(neutro_20mix_ME~neutro_20mix_mod[,g], main="20% neutro")
}

#some genes separate convolutions well...others don't

# 2...
# repeat loop for each mixture

neutro90_ME_dist = matrix(0, nrow = ncol(neutro_90mix_mod), ncol = ncol(neutro_90mix_mod))
for(i in 1:ncol(neutro_90mix_mod)){
  gene_ME_dist = matrix(0, nrow = ncol(neutro_90mix_mod), ncol = nrow(neutro_90mix_mod))
  for(j in 1:nrow(neutro_90mix_mod)){
    gene_scaler = abs(neutro_90mix_mod[j,i]-neutro_90mix_ME[j])
    ME_dist = sapply(unlist(neutro_90mix_mod[j,]), function(x) abs(x - neutro_90mix_ME[j])/gene_scaler)
    gene_ME_dist[,j] = ME_dist
  }
  neutro_avg_dist = apply(gene_ME_dist, 1, mean)
  neutro90_ME_dist[i,] <- neutro_avg_dist
}

# Flipping the above loop to get sample-specific distances per gene
neutro10_ME_dist = array(matrix(0, nrow = ncol(neutro_10mix_mod), ncol = ncol(neutro_10mix_mod)),
                          dim = c(ncol(neutro_10mix_mod),ncol(neutro_10mix_mod),nrow(neutro_10mix_mod)))
for(i in 1:nrow(neutro_10mix_mod)){
  # browser()
  gene_ME_dist = matrix(0, nrow = ncol(neutro_10mix_mod), ncol = ncol(neutro_10mix_mod))
  for(j in 1:ncol(neutro_10mix_mod)){
    # gene_scaler = neutro_10mix_mod[i,j]-neutro_10mix_ME[i]
    ME_dist = sapply(unlist(neutro_10mix_mod[i,]), function(x) x/neutro_10mix_mod[i,j])
    gene_ME_dist[j,] = ME_dist
  }
  # neutro_avg_dist = apply(gene_ME_dist, 1, mean)
  neutro10_ME_dist[,,i] <- gene_ME_dist
}

par(mfrow=c(6,3))
par(mar=c(1,1,1,1))
neutroMixList = list('training' = neutro_ME_dist, '100' = neutro100_ME_dist, '80' = neutro80_ME_dist,
                     '60' = neutro60_ME_dist, '40' = neutro40_ME_dist, '20' = neutro20_ME_dist)

randomGenes <- sample(1:sum(neutro_colors=='orange'),3)
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

# Trying to look at genes by individual samples
# randomGenes <- sample(1:sum(neutro_colors=='darkturquoise'),3)
neutroArrayList = list('100' = neutro100_ME_dist, '80' = neutro80_ME_dist,
                     '60' = neutro60_ME_dist, '40' = neutro40_ME_dist, '20' = neutro20_ME_dist)
randomSample <- sample(1:100,1)
randomGenes <- sample(1:sum(neutro_colors=='orange'),3)
for(g in randomGenes){
  x = neutro_ME_dist[g,]
  x.fit <- seq(min(x), max(x), length=length(x))
  h = hist(neutro_ME_dist[g,], breaks=20, xlim = c(0.5,2), main = 'training')
  fit = fitdistr(neutro_ME_dist[g,], 'Gamma')
  x.gamma <- dgamma(x.fit, shape=fit$estimate[1], rate=fit$estimate[2])
  lines(x.fit, x.gamma/max(x.gamma)*max(h$density), col='green', lwd=2)
}
for(data in 1:length(neutroArrayList)){
  for(g in randomGenes){
    x = unname(neutroArrayList[[data]])[g,,randomSample]
    x.fit <- seq(min(x), max(x), length=length(x))
    h = hist(neutroArrayList[[data]][g,,randomSample], breaks=20, xlim = c(0.5,2), main = names(neutroArrayList[data]))
    fit = fitdistr(neutroArrayList[[data]][g,,randomSample], 'Gamma')
    x.gamma <- dgamma(x.fit, shape=fit$estimate[1], rate=fit$estimate[2])
    lines(x.fit, x.gamma/max(x.gamma)*max(h$density), col='green', lwd=2)
  }
}


#####
# 1a. Calculate a divergence/entropy score from fitted distributions
# 1b. Divergence doesn't work so well, a correlation/covariance seems to match much better
#     ...interestingly, the gene-wise covariance don't seem to be that tightly correlated
# 2. Make a full correlation data frame of all mixtures
#     - Make a train and test set
# 3. Model the covariance with mixture proporation
#####

randomGene <- sample(1:sum(neutro_colors=='orange'),1)
randomSample <- sample(1:100,1)
for(data in 1:length(neutroArrayList)){
  # gene_divergence = sapply(1:sum(neutro_colors=='darkturquoise'),
                          # function(x) KL.divergence(neutro_ME_dist[,x],neutroArrayList[[data]][,x,18],2)[2])
  # print(mean(gene_divergence))
  # print(KL.divergence(neutro_ME_dist[randomGene,],neutroArrayList[[data]][randomGene,,randomSample]))
  # print(cor(neutro_ME_dist[randomGene,],neutroArrayList[[data]][randomGene,,randomSample]))
  gene_correlation = sapply(1:sum(neutro_colors=='orange'),
                            function(x) cor(neutro_ME_dist[x,],neutroArrayList[[data]][x,,randomSample]))
  print(mean(gene_correlation))
  # plot(neutro_ME_dist[,randomGene],neutroArrayList[[data]][,randomGene,25])
}

cov_test =t(sapply(1:3, function(i) sapply(1:5,
       function(x) cov(neutro_ME_dist[x,],neutroArrayList[[5]][x,,i]))))


# Step 2...
neutroArrayList = list('100' = neutro100_ME_dist, '90' = neutro90_ME_dist, '80' = neutro80_ME_dist,
                       '70' = neutro70_ME_dist, '60' = neutro60_ME_dist, '50' = neutro50_ME_dist,
                       '40' = neutro40_ME_dist, '30' = neutro30_ME_dist, '20' = neutro20_ME_dist,
                       '10' = neutro10_ME_dist, '0' = neutro0_ME_dist)

initial = 0
for(idx in 1:length(neutroArrayList)){
  if(initial == 0){
    all_cov_matrix = t(sapply(1:100, function(i) sapply(1:sum(neutro_colors=='orange'), 
                                                  function(x) cov(neutro_ME_dist[x,],neutroArrayList[[idx]][x,,i]))))
    all_cov_matrix = cbind(all_cov_matrix, rep(as.integer(names(neutroArrayList)[idx]),100))
    colnames(all_cov_matrix) = c(colnames(neutro_orange_df), "mix")
    initial = 1
  } else {
    cov_matrix = t(sapply(1:100, function(i) sapply(1:sum(neutro_colors=='orange'), 
                                                        function(x) cov(neutro_ME_dist[x,],neutroArrayList[[idx]][x,,i]))))
    cov_matrix = cbind(cov_matrix, rep(as.integer(names(neutroArrayList)[idx]),100))
    colnames(cov_matrix) = c(colnames(neutro_orange_df), "mix")
    all_cov_matrix = rbind(all_cov_matrix, cov_matrix)
    }
  }


test_idx <- sample(1:1000,100)
test_cov_matrix <- all_cov_matrix[test_idx,]
train_cov_matrix <- all_cov_matrix[-test_idx,]

# Step 3...

#Finding the best alpha

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


# somewhere around the middle is good, but alpha = 0.9 was best
# honestly, 0.1-0.9 are all pretty okay
# It seems more accurate for the lower percentages
# for now I'll try 3...

#NEED TO INCLUDE 0% SAMPLES TOO!!!
# testing a 0%

neutro_0mix = read.csv('../mixed_dfs/neutro_coarse_all_0.csv')
# a little restructring...
neutro_0mix <- as.data.frame(neutro_0mix)
row.names(neutro_0mix) = neutro_0mix$gene_symbol_sql
neutro_0mix = subset(neutro_0mix, select = -c(gene_symbol_sql))
neutro_0mix = t(neutro_0mix)
# a little more formatting
mix_cols = colnames(neutro_0mix)[colnames(neutro_0mix) %in% colnames(coarse_df)]
neutro_0mix = neutro_0mix[,mix_cols]
missing_genes = colnames(coarse_df)[(!(colnames(coarse_df) %in% colnames(neutro_0mix)))]
median_expression = apply(neutro_0mix,1,median)
missing_genes_df = sapply(missing_genes, function(g) median_expression, simplify = 'array')
neutro_0mix_all = cbind(neutro_0mix, missing_genes_df)
neutro_0mix_all = neutro_0mix_all[,colnames(coarse_df)]
neutro_0mix_mod = neutro_0mix_all[,neutro_colors=='orange']

neutro0_ME_dist = array(matrix(0, nrow = ncol(neutro_0mix_mod), ncol = ncol(neutro_0mix_mod)),
                         dim = c(ncol(neutro_0mix_mod),ncol(neutro_0mix_mod),nrow(neutro_0mix_mod)))
for(i in 1:nrow(neutro_0mix_mod)){
  gene_ME_dist = matrix(0, nrow = ncol(neutro_0mix_mod), ncol = ncol(neutro_0mix_mod))
  for(j in 1:ncol(neutro_0mix_mod)){
    ME_dist = sapply(unlist(neutro_0mix_mod[i,]), function(x) x/neutro_0mix_mod[i,j])
    gene_ME_dist[j,] = ME_dist
  }
  neutro0_ME_dist[,,i] <- gene_ME_dist
}

mix0_cov_matrix = t(sapply(1:100, function(i) sapply(1:sum(neutro_colors=='orange'), 
                                                    function(x) cov(neutro_ME_dist[x,],neutro0_ME_dist[x,,i]))))
mix0_cov_matrix = cbind(mix0_cov_matrix, rep(0,100))
colnames(mix0_cov_matrix) = c(colnames(neutro_orange_df), "mix")

yhat0 <- predict(fit3, s=fit3$lambda.1se, newx=subset(mix0_cov_matrix, select = -c(mix)))


# for predicting a single sample
# predict(mono_ME_model, s=mono_ME_model$lambda.1se, newx = t(as.matrix(unlist(mono_steelblue_df[3,]))))

save(b_avg_dist, cd4_avg_dist, cd8_avg_dist, endo_avg_dist, fibro_avg_dist, mono_avg_dist, neutro_avg_dist, nk_avg_dist,
     b_mod_geneNames, cd4_mod_geneNames, cd8_mod_geneNames, endo_mod_geneNames, fibro_mod_geneNames, mono_mod_geneNames,
     neutro_mod_geneNames, nk_mod_geneNames, b_ME_model, cd4_ME_model, cd8_ME_model, endo_ME_model, fibro_ME_model,
     mono_ME_model, neutro_ME_model, nk_ME_model, file='WGCNA_module_output.RData')

