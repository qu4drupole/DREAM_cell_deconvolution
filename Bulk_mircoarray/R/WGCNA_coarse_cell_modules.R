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
# testing 2. silhouette score with neutro vs. b cells
#####
b_cell_df = coarse_df[cellTypes=='B.cells',]

#euclidean norm of 'grey' module
# calculating the Frobenius norm which is essentially the euclidean norm...?
neutro_mod_test = t(neutro_cell_df[,neutro_colors=='grey'])
b_mod_test = t(b_cell_df[,neutro_colors=='grey'])

norm(neutro_mod_test, "F") #627.96
norm(b_mod_test, "F") #1129.926



binomialCoef <- function(n,k){
  #these numbers are too large to use factorial
  value <- lgamma(n+1)-((lgamma(k+1)+lgamma(n-k+1)))
  value <- exp(value)
  return(value)
}

# 1...
all_ct_modules = list('b' = b_colors, 'cd4' = cd4_colors, 'cd8' = cd8_colors, 'nk' = nk_colors,
                      'neutro' = neutro_colors, 'mono' = mono_colors, 'fibro' = fibro_colors, 'endo' = endo_colors)

ct_modules = unlist(all_ct_modules['cd8'])
test_modules = all_ct_modules[-c(names(all_ct_modules)=='cd8')]
ct_module_probs = list()

# 2...
for(color in unique(ct_modules)){
  module_total_p = 0
  module_genes = rep(0,ncol(coarse_df))
  module_genes[ct_modules==color] = 1
  
  if(sum(module_genes) > 500){
    next
  }
  
  for(ct in names(test_modules)){
    for(mod in unique(unlist(test_modules[ct]))){
      # browser()
      test_genes = rep(0,ncol(coarse_df))
      test_genes[unlist(test_modules[ct])==mod] = 1
      
      K = sum((module_genes+test_genes)==2)
      k = floor(K/2)
      N = sum(test_genes == 1)
      n = floor(N/2)
      
      if(K == 0 | N > 500){ 
        p = 0
      } else {
        p = (binomialCoef(K,k)*binomialCoef(N-K,n-k))/binomialCoef(N,n)  
      }
      module_total_p = module_total_p + p
    }
  }
  
  ct_module_probs[paste(color)] = module_total_p
}

# ct_module_probs[order(unlist(ct_module_probs))][1:5]
# sum(ct_modules=='yellow')

plot(b_cyan_df[,10],b_MEs$MEcyan)

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

