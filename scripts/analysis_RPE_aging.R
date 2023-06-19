##########################################################################
##########################################################################
# Project: RPE cell line analysis for Birgit 
# Script purpose: processing published microarray data and RNA-seq data to identify Birgit's cell line's age 
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue May  2 15:54:33 2023
##########################################################################
##########################################################################
rm(list = ls())

version.analysis = '_RPE_aging_20230505'

resDir = paste0("../results/aging", version.analysis)
RdataDir = paste0('../results/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

## packages and functions
require(ggplot2)
require("sva")
require(limma)
require(preprocessCore)

make.pca.plots = function(fpm, ntop = 1000, pca.dim = c(1, 2), scale.data.pca = FALSE)
{
  library(factoextra)
  
  xx = as.matrix(fpm)
  vars = apply(xx, 1, var)
  xx = xx[order(-vars), ]
  xx = xx[1:min(ntop, nrow(xx)), ]
  #par(mfrow = c(2,2))
  #pairs(xx[, c(1:4)])
  #plot.pair.comparison.plot(xx[, c(1:4)], linear.scale = FALSE)
  #plot.pair.comparison.plot(xx[, c(9:16)], linear.scale = FALSE)
  
  res.pca <- prcomp(t(xx), scale = scale.data.pca)
  #res.var <- get_pca_var(res.pca)
  
  fviz_pca_ind(res.pca,
               axes = pca.dim,
               col.ind = "cos2", # Color by the quality of representation
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE     # Avoid text overlapping
  )
  
}

########################################################
########################################################
# Section I: process RNA-seq data from Bulter_2021_RNAseq_GSE159435
# 
########################################################
########################################################
xlist = list.files(path =  '../Bulter_2021_RNAseq_GSE159435_RAW', pattern = '*.txt', full.names = TRUE)

for(n in 1:length(xlist))
{
  # n = 1
  test = read.table(file = xlist[n])
  if(n == 1){
    counts = test
  }else{
    
    counts = cbind(counts, test[match(counts[,1], test[,1]), 2])
    
  }
}

colnames(counts) = c('geneID', basename(xlist))
colnames(counts) = gsub('_filtered_count.txt', '', colnames(counts))

# metadata from https://www.ncbi.nlm.nih.gov/geo/geo2r/?acc=GSE159435
design = data.frame(GEO.number =  sapply(colnames(counts)[-1], 
                                         function(x) unlist(strsplit(as.character(x), '_'))[1]),
                    donor = sapply(colnames(counts)[-1], 
                                   function(x) unlist(strsplit(as.character(x), '_'))[2]))
design$batch = c(rep('s1', 7), rep('s2', 6))
design$age = c(31, 40, 51, 67, 77, 86, 93, 61, 62, 71, 73, 73, 92)
design$sex = c('M', 'F', 'M', 'F', 'M', 'M', 'F', 'F', 'M', 'M', 'M', 'F', 'M')

annots = read.delim('../data/Human.GRCh38.p13.annot.tsv')
counts = counts[, c(2:14, 1)]

mm = match(counts$geneID, annots$EnsemblGeneID)
counts = data.frame(counts, annots[mm, c(2, 1, 4, 5, 6)], stringsAsFactors = FALSE)

counts = counts[which(!is.na(counts$Symbol)), ]
rownames(counts) = counts$Symbol

counts = counts[, c(1:13)]

save(counts, design, file = paste0(RdataDir, 'design_countTable_RPE_RNAseq_Bulter_2021.Rdata'))

## DESeq2 normalization and PCA
require(DESeq2)


dds <- DESeqDataSetFromMatrix(as.matrix(counts), DataFrame(design), design = ~ 1)

ss = rowSums(counts(dds))

hist(log10(ss), breaks = 100);abline(v = log10(20), lwd = 2.0, col = 'red')
cat(length(which(ss>50)), ' gene selected \n')

dds = dds[which(ss>50), ]

dds = estimateSizeFactors(dds)

vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

pca=plotPCA(vsd, intgroup = c('age', 'batch', 'sex'), returnData = FALSE)
print(pca)

pca2save = as.data.frame(plotPCA(vsd, intgroup = c('age', 'batch', 'sex', 'donor'), 
                                 returnData = TRUE, ntop = 500))
pca2save$name = paste0(pca2save$donor, '_', pca2save$age, '_', pca2save$sex)
pca2save$age = factor(pca2save$age)

ggplot(data=pca2save, aes(PC1, PC2, label = name, color= age, shape = batch))  + 
  geom_point(size=3) + 
  geom_text(hjust = 1.2, nudge_y = 1.5, size=2.5)

ggsave(paste0(resDir, '/RPE_RNAseq.pdf'), 
       width=12, height = 8)

save(dds, design, file = paste0(RdataDir, 'dds_design_RPE_RNAseq_Bulter_2021.Rdata'))

load(file = paste0(RdataDir, 'dds_design_RPE_RNAseq_Bulter_2021.Rdata'))
cpm = fpm(dds)
cpm = log2(cpm + 2^-4)

xx = apply(cpm, 1, function(x){
  t = as.numeric(design$age)
  fit = lm(x ~ t)
  beta = summary(fit)$coefficients[2,1]
  pval = summary(fit)$coefficients[2,4]
  intercept = summary(fit)$coefficients[1,1]
  return(c(beta, pval, intercept))
})


xx = t(xx)
xx = xx[, c(2,1,3)]
colnames(xx) = c('pval', 'beta', 'intercept') # slope and offset

cpm = data.frame(cpm, xx, stringsAsFactors = FALSE)
cpm = cpm[order(cpm$pval), ]

save(dds, design, cpm, file = paste0(RdataDir, 'dds_design_ageGenes_RPE_RNAseq_Bulter_2021.Rdata'))

load(file = paste0(RdataDir, 'dds_design_ageGenes_RPE_RNAseq_Bulter_2021.Rdata'))

cpm[grep('LRAT', rownames(cpm)),]

pdf(paste0(resDir, "/RNAseq_geneExp_age_lm_fitting.pdf"), height = 8, width =10)

for(n in 1:500)
{
  cat('n -- ', n, '\n')
  x =  as.numeric(cpm[n, c(1:13)])
  t =  as.numeric(design$age)
  fit = lm(x ~ t)
  
  plot(t, x, cex = 1.2, col = 'darkblue', main = paste0(rownames(cpm)[n],  ' : pval -- ', 
                                                        signif(cpm$pval[n], d =2)))
  abline(fit$coefficients, lwd = 1.2, col = 'red')
  
}

dev.off()


########################################################
########################################################
# Section II : process micorarray data from GSE18811
# 
########################################################
########################################################

# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0, DESeq2 1.38.3
################################################################
#   Data plots for selected GEO samples
library(GEOquery)
library(limma)
library(umap)

# load series and platform data from GEO 
gset <- getGEO("GSE18811", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# microarray expression matrix
ex <- exprs(gset)
#ex[which(ex <= 0)] <- NaN

# sample info and select the relevant samples after dicussion with Birgit
design = gset@phenoData@data[, c(10, 38,39)]

kk = grep('native fetal RPE|adult native RPE', design$characteristics_ch1)
design = design[kk, ]
ex = ex[,kk]
design$group = NA
design$group[grep('fetal', design$`cell type:ch1`)] = 'fetal'
design$group[grep('adult', design$`cell type:ch1`)] = 'adult'
design$age = design$`gestational age:ch1`
design$age = gsub(' weeks', 'w', design$age)
design$age[c(4:6,10)] = c('76', '78', '79', '78')
design$condition = paste0(design$group, '.', design$age)
colnames(ex) = paste0(design$age, '_', colnames(ex))

## prob id to gene mapping
#require('hgu133plus2.db')
#fdata(gset)
#colnames(fData(data))
library(openxlsx)
#mapping = read.csv('../data/GPL570-55999.csv', sep = '\t')
mapping = read.xlsx('../data/GPL570-55999_v1.xlsx')
mapping = mapping[, c(1, 2, 10, 11, 12)]
mapping = mapping[, c(1, 4)]

mm = match(rownames(ex), mapping$ID)
length(which(!is.na(mm)))
mapping = mapping[mm, ]
mapping$Gene.Symbol = gsub(' /// ', '_', mapping$Gene.Symbol)

save(mapping, design, ex, file = paste0(RdataDir, 'design_probAnnot_expr.Rdata'))

##########################################
# processing the microarray data
##########################################
mat = as.matrix(ex)

mat = log2(mat)

## quantil normalization
library(preprocessCore)
mat.norm = normalize.quantiles(mat, copy = TRUE, keep.names = TRUE)

#mat.norm = limma::normalizeBetweenArrays(mat, )
colnames(mat.norm) = colnames(mat)

## construct gene-to-sample matrix
ggs = unique(mapping$Gene.Symbol)
ggs = ggs[which(!is.na(ggs))]

res = matrix(NA, ncol = ncol(mat.norm), nrow = length(ggs))
colnames(res) = colnames(mat.norm)
rownames(res) = ggs

for(n in 1:nrow(res))
{
  if(n%%500 == 0) cat(n, '\n')
  jj = which(mapping$Gene.Symbol == rownames(res)[n])
  if(length(jj) > 1) {
    res[n, ] = apply(mat.norm[jj, ], 2, median)
  }else{
    res[n, ] = mat.norm[jj, ]
  }
}

f <- factor(design$group, levels=c("fetal","adult"))
cc <- model.matrix(~0+f)
colnames(cc) <- c("fetal","adult")

fit <- lmFit(res, cc)
contrast.matrix <- makeContrasts(adult - fetal, levels=cc)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

#A list of top genes for RNA2 versus RNA1 can be obtained from
topTable(fit2, coef=1, adjust="BH")

tops = topTable(fit2, adjust="BH", number = nrow(res),  genelist = rownames(res))[, c(2, 5,6)]
colnames(tops) = paste0(colnames(tops), '_adult.vs.fetal')
res = data.frame(res, tops[match(rownames(res), rownames(tops)), ])

#The outcome of each hypothesis test can be assigned using
#results <- decideTests(fit2)
#xx = data.frame(fit2$fit2$p.value)
#colnames(xx) = c('adult.vs.fetal')
#xx$pval.max = apply(as.matrix(-log10(xx)), 1, max)

#res = data.frame(res, xx)
res = res[order(res$P.Value_adult.vs.fetal), ]

save(mapping, design, ex, res, 
     file = paste0(RdataDir, 'design_probAnnot_expr_gene.sample.matrix.Rdata'))

##########################################
# compare RNA-seq age-genes with microarray RPE
##########################################
load(file = paste0(RdataDir, 'dds_design_ageGenes_RPE_RNAseq_Bulter_2021.Rdata'))
load(file = paste0(RdataDir,  'Bulter_2021_RPE_ageGenes_pval0.01_prediction.allTimepoints_design.Rdata'))
genes2keep = readRDS(file = paste0(RdataDir, 'genes_to_keep_RNAseq_vs_RPEChoroid.microarray.rds'))

kk = which(cpm$pval<0.01); cat(length(kk), ' age-related genes \n')
xx = cpm[kk, c(1:13)]
colnames(xx) = paste0(design$donor, '_', design$sex, '_', design$age)
rm(design)

xx0 = xx0[match(genes2keep, rownames(xx0)), ]
xx = xx[match(genes2keep, rownames(xx)), ]

load(file = paste0(RdataDir, 'design_probAnnot_expr_gene.sample.matrix.Rdata'))


pdf(paste0(resDir, "/RNAseq_ageGenes_vs_fetal.Old.microarray.pdf"), height = 8, width =10)

#jj_macular = grep('_macular_', colnames(res))
tt_microarray = design$age
tt_microarray[grep('w', tt_microarray)] = 0.30
tt_microarray = as.numeric(tt_microarray)
#jj_extra = grep('_extramacular_', colnames(res))
#tt_extra = as.numeric(design$age[jj_extra])

slopes = c()

for(n in 1:nrow(xx0))
{
  # n = 1
  ii = which(rownames(res) == rownames(xx0)[n])
  test = c(cpm$beta[which(rownames(cpm) == rownames(xx0)[n])])
  
  if(length(ii) == 1){
    cat(n, '-- gene : ', rownames(xx0)[n],   '--\n')
    x =  as.numeric(xx0[n, -c(1:13)])
    x = as.numeric(scale(x, center = TRUE, scale = FALSE))
    t =  as.numeric(design0$age[-c(1:13)])
    x1 = as.numeric(scale(as.numeric(res[ii, 1:10]), center = TRUE, scale = FALSE))
    #x2 = as.numeric(scale(res[ii, jj_extra], center = TRUE, scale = FALSE))
    plot(t, x, cex = 1.0, type = 'p', ylim = range(c(x, x1, x2)),
         main = paste0(rownames(xx0)[n]), xlab = 'age (year)', 
         ylab = 'centered log2(Expression)', col = 'red')
    points(t, x, lwd = 1.2, col = 'red', type = 'l')
    #x =  as.numeric(xx0[n, c(1:13)])
    #t =  as.numeric(design0$age[c(1:13)])
    #points(t, x, lwd = 1.2, col = 'blue', type = 'p')
    
    points(tt_microarray, x1, cex = 1.2, col = 'darkorange')
    #points(tt_extra, x2, cex = 1.2, col = 'black')
    fit1 = lm(x1 ~ tt_microarray)
    test = c(test, coefficients(fit1)[2])
    #abline(coefficients(fit1), lwd = 1.2, col = 'darkorange')
    #fit2 = lm(x2 ~ tt_extra)
    #abline(coefficients(fit2), lwd = 1.2, col = 'black')
    slopes = rbind(slopes, test)
  }else{
    cat(n, '-- gene : ', rownames(xx0)[n],   'not found--\n')
    slopes = rbind(slopes, c(test, NA))
  }
  
}

dev.off()

sels = which(!is.na(slopes[,2]) & sign(slopes[,1]) == sign(slopes[,2]))
genes2keep = rownames(xx0)[sels]

saveRDS(genes2keep, 
        file = paste0(RdataDir, 'genes_to_keep_RNAseq_vs_RPEChoroid.microarray_RPE.fetal.age.microarray.rds'))

########################################################
########################################################
# Section II.5:
# download and process microarray data from 
# Newman AM, Gallo NB, Hancox LS, Miller NJ et al. 
# Genome Med 2012 Feb 24;4(2):16. 
########################################################
########################################################
library(GEOquery)
library(limma)
library(umap)

# load series and platform data from GEO
gset <- getGEO("GSE29801", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL4133", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# microarray expression matrix
ex <- exprs(gset)

# sample info and select the relevant samples after dicussion with Birgit
design = gset@phenoData@data
kk = grep('RPE-choroid', design$`tissue:ch1`)
design = design[kk, ]
ex = ex[,kk]

kk = which(design$`ocular disease:ch1` == 'normal')
design = design[kk, ]
ex = ex[,kk]

design = data.frame(design$title, design$geo_accession, design$platform_id, 
                    design$`age (years):ch1`,
                    #design$`ocular disease:ch1`,
                    design$`patient identifier:ch1`, design$`tissue:ch1`)
colnames(design)[1:6] = c('sample', 'accession', 'platform.id', 'age', 'patient.id', 'tissue')
design$tissue = gsub(' RPE-choroid', '', design$tissue)

design$group = paste0(design$accession, '_', design$tissue, '_', design$age)

colnames(ex) = design$group

## prob id to gene mapping
mapping = read.csv2('../Newman_2012/GSE29801_RAW/GPL4133_old_annotations.csv')
mm = match(rownames(ex), mapping$ID)
length(which(is.na(mm)))

mapping = mapping[mm, ]

#mapping = read.xlsx('../data/GPL570-55999_v1.xlsx')
mapping = mapping[, c(1, 4, 6, 10)]
#mapping = mapping[, c(1, 4)]

#mm = which(mapping$GENE_SYMBOL != '')
#length(which(!is.na(mm)))
#mapping = mapping[mm, ]
#mapping$Gene.Symbol = gsub(' /// ', '_', mapping$Gene.Symbol)

save(mapping, design, ex, file = paste0(RdataDir, 'Newman_2012_design_probAnnot_expr.Rdata'))

##########################################
# processing the microarray data
##########################################
mat = as.matrix(ex)

mat = log2(mat)

## quantil normalization
library(preprocessCore)
mat.norm = normalize.quantiles(mat, copy = TRUE, keep.names = TRUE)

#mat.norm = limma::normalizeBetweenArrays(mat, )
colnames(mat.norm) = colnames(mat)

## construct gene-to-sample matrix
ggs = unique(mapping$GENE_SYMBOL)
ggs = ggs[which(!is.na(ggs) & ggs != '')]

res = matrix(NA, ncol = ncol(mat.norm), nrow = length(ggs))
colnames(res) = colnames(mat.norm)
rownames(res) = ggs

for(n in 1:nrow(res))
{
  if(n%%500 == 0) cat(n, '\n')
  jj = which(mapping$GENE_SYMBOL == rownames(res)[n])
  
  if(length(jj) > 1) {
    res[n, ] = apply(mat.norm[jj, ], 2, median)
  }else{
    res[n, ] = mat.norm[jj, ]
  }
}


save(design, res, 
     file = paste0(RdataDir, 'Newman_et_al_2012_design_normExpr.Rdata'))

make.pca.plots(res[, grep('_extramacular_', colnames(res))], ntop = 1000, pca.dim = c(1, 3))

##########################################
# compare the RNA-seq and micro-array data in Newman 2012
##########################################
load(file = paste0(RdataDir, 'dds_design_ageGenes_RPE_RNAseq_Bulter_2021.Rdata'))
load(file = paste0(RdataDir,  'Bulter_2021_RPE_ageGenes_pval0.01_prediction.allTimepoints_design.Rdata'))

kk = which(cpm$pval<0.01); cat(length(kk), ' age-related genes \n')
xx = cpm[kk, c(1:13)]
colnames(xx) = paste0(design$donor, '_', design$sex, '_', design$age)
rm(design)

load(file = paste0(RdataDir, 'Newman_et_al_2012_design_normExpr.Rdata'))

## compare the microarray and RNA-seq
Compare_RNAseq_micorray = FALSE
if(Compare_RNAseq_micorray){
  library(preprocessCore)
  yy1 = cpm[, c(1:13)]
  yy2 = res
  
  ggs_intersect = intersect(rownames(yy1), rownames(yy2))
  
  yy1 = yy1[match(ggs_intersect, rownames(yy1)), ]
  yy2 = yy2[match(ggs_intersect, rownames(yy2)), ]
  
  yy1 = apply(yy1, 1, median)
  #yy2 = apply(yy2, 1, median)
  
  yy3 = preprocessCore::normalize.quantiles.use.target(as.matrix(yy2), target = as.vector(yy1))
  colnames(yy3) = colnames(yy2)
  rownames(yy3) = rownames(yy2)
  res = yy3
  
  #yy3 = apply(yy3, 1, median)
  
  plot(yy1, yy3, cex = 0.5)
  abline(0, 1, lwd = 1.5, col = 'red')
  
  
}

pdf(paste0(resDir, "/RNAseq_ageGenes_pval0.01_vs_Newman.2012_quantileNormalization.pdf"), 
    height = 8, width =10)

jj_macular = grep('_macular_', colnames(res))
tt_macular = as.numeric(design$age[jj_macular])
jj_extra = grep('_extramacular_', colnames(res))
tt_extra = as.numeric(design$age[jj_extra])


slopes = c()
for(n in 1:nrow(xx0))
{
  # n = 1
  test = c(cpm$beta[which(rownames(cpm) == rownames(xx0)[n])])
  
  ii = which(rownames(res) == rownames(xx0)[n])
  if(length(ii) == 1){
    cat(n, '-- gene : ', rownames(xx0)[n],   '--\n')
    x =  as.numeric(xx0[n, -c(1:13)])
    x = as.numeric(scale(x, center = TRUE, scale = FALSE))
    t =  as.numeric(design0$age[-c(1:13)])
    x1 = as.numeric(scale(as.numeric(res[ii, jj_macular]), center = TRUE, scale = FALSE))
    x2 = as.numeric(scale(as.numeric(res[ii, jj_extra]), center = TRUE, scale = FALSE))
    plot(t, x, cex = 1.0, type = 'p', ylim = range(c(x, x1, x2)),
         main = paste0(rownames(xx0)[n]), xlab = 'age (year)', 
         ylab = 'centered log2(Expression)', col = 'red')
    points(t, x, lwd = 1.2, col = 'red', type = 'l')
    #x =  as.numeric(xx0[n, c(1:13)])
    #t =  as.numeric(design0$age[c(1:13)])
    #points(t, x, lwd = 1.2, col = 'blue', type = 'p')
   
    points(tt_macular, x1, cex = 1.2, col = 'darkorange')
    points(tt_extra, x2, cex = 1.2, col = 'black')
    fit1 = lm(x1 ~ tt_macular)
    abline(coefficients(fit1), lwd = 1.2, col = 'darkorange')
    fit2 = lm(x2 ~ tt_extra)
    test = c(test, coefficients(fit1)[2], coefficients(fit2)[2])
    abline(coefficients(fit2), lwd = 1.2, col = 'black')
    slopes = rbind(slopes, test)
    
  }else{
    cat(n, '-- gene : ', rownames(xx0)[n],   'not found--\n')
    slopes = rbind(slopes, c(test, NA, NA))
  }
  
}

dev.off()

## filter the age-dependent genes
rownames(slopes) = rownames(xx0)
sels = which(!is.na(slopes[,2]) & !is.na(slopes[,3]))
slopes = slopes[sels, ]
sels = c()
sels = apply(slopes, 1, function(x) {sign(x[1]) == sign(x[2]) & sign(x[1]) == sign(x[3])})
genes2keep = rownames(slopes)[sels]

saveRDS(genes2keep, file = paste0(RdataDir, 'genes_to_keep_RNAseq_vs_RPEChoroid.microarray.rds'))

yy = data.frame(xx[match(genes2keep, rownames(xx)), ], 
                res[match(genes2keep, rownames(res)), ])

pca2save <- prcomp(t(as.matrix(yy)), scale = TRUE)
#pca2save = pca2save[, -c(1:2)]
pca2save = data.frame(pca2save$x, design0)
pca2save = data.frame(pca2save$x)
pca2save$name = rownames(pca2save)
pca2save$age = factor(pca2save$age)

library(ggplot2)
ggplot(data=pca2save, aes(PC1, PC2, label = name))  + 
  geom_point(size=3) + 
  geom_text(hjust = 0.3, nudge_y = 0.3, size=2.5)

#ggs2keep = readRDS(file = paste0(RdataDir, 'genes_to_keep_RNAseq_vs_RPEChoroid.microarray.rds'))

##########################################
# select a list of RPE aging-related genes
##########################################
load(file = paste0(RdataDir, 'dds_design_ageGenes_RPE_RNAseq_Bulter_2021.Rdata'))

kk = which(cpm$pval<0.01); cat(length(kk), ' age-related genes \n')

xx = cpm[kk, c(1:13)]
colnames(xx) = paste0(design$donor, '_', design$sex, '_', design$age)

#make.pca.plots(xx, pca.dim = c(1, 2), ntop = 200)
pca2save <- prcomp(t(as.matrix(xx)), scale = TRUE)
#pca2save = pca2save[, -c(1:2)]
pca2save = data.frame(pca2save$x, design)

pca2save$name = paste0(pca2save$donor, '_', pca2save$age, '_', pca2save$sex)
pca2save$age = factor(pca2save$age)

ggplot(data=pca2save, aes(PC1, PC2, label = name, color= age, shape = batch))  + 
  geom_point(size=3) + 
  geom_text(hjust = 0.3, nudge_y = 0.3, size=2.5)

ggsave(paste0(resDir, '/RPE_RNAseq_PCA_910age.genes.pdf'), 
       width=8, height = 4)


## add extrapolated data
kk = which(cpm$pval<0.01); cat(length(kk), ' age-related genes \n')
xx = cpm[kk, c(1:13)]
colnames(xx) = paste0(design$donor, '_', design$sex, '_', design$age)

tt = seq(0, 90, by = 10)
xx0 = matrix(NA, nrow = nrow(xx), ncol = length(tt))
rownames(xx0) = rownames(xx)
colnames(xx0) = paste0('syn_age_', tt)

for(n in 1:nrow(xx0))
{
  cat(n, ' --- gene : ', rownames(xx0)[n], '--\n')
  a = cpm$intercept[which(rownames(cpm) == rownames(xx0)[n])]
  b = cpm$beta[which(rownames(cpm) == rownames(xx0)[n])]
  xx0[n, ] = a + b*tt
}

xx = cbind(xx, xx0)
design0 = design[, c(3:5)]
design0 = rbind(design0, data.frame(batch = rep('synthetic', length(tt)), 
                                    age = tt, 
                                    sex = "M.F"))
rownames(design0) = colnames(xx)
xx0 = xx
save(xx0, design0,  file = paste0(RdataDir, 
                                  'Bulter_2021_RPE_ageGenes_pval0.01_prediction.allTimepoints_design.Rdata'))

pca2save <- prcomp(t(as.matrix(xx)), scale = TRUE)
#pca2save = pca2save[, -c(1:2)]
pca2save = data.frame(pca2save$x, design0)

pca2save$name = rownames(pca2save)
pca2save$age = factor(pca2save$age)

ggplot(data=pca2save, aes(PC1, PC2, label = name, color= age, shape = batch))  + 
  geom_point(size=3) + 
  geom_text(hjust = 0.3, nudge_y = 0.3, size=2.5)

ggsave(paste0(resDir, '/RPE_RNAseq_PCA_age.genes_pval0.01_predition_alltimepoints.pdf'), 
       width=8, height = 4)


########################################################
########################################################
# Section IV: 
# test to integrate RNA-seq data and microarray 
# 
########################################################
########################################################
load(file = paste0(RdataDir, 'dds_design_ageGenes_RPE_RNAseq_Bulter_2021.Rdata'))
load(file = paste0(RdataDir,  'Bulter_2021_RPE_ageGenes_pval0.01_prediction.allTimepoints_design.Rdata'))

kk = which(cpm$pval<0.01); cat(length(kk), ' age-related genes \n')
xx = cpm[kk, c(1:13)]
colnames(xx) = paste0(design$donor, '_', design$sex, '_', design$age)
rm(design)

#load(file = paste0(RdataDir, 'Newman_et_al_2012_design_normExpr.Rdata'))
load(file = paste0(RdataDir, 'design_probAnnot_expr_gene.sample.matrix.Rdata'))
#yy = data.frame(yy, res[match(genes2keep, rownames(res)), c(1:10)])

genes2keep = readRDS(file = paste0(RdataDir, 
                                   'genes_to_keep_RNAseq_vs_RPEChoroid.microarray_RPE.fetal.age.microarray.rds'))


## compare the microarray and RNA-seq
Compare_RNAseq_micorray = FALSE
if(Compare_RNAseq_micorray){
  library(preprocessCore)
  yy1 = cpm[, c(1:13)]
  colnames(yy1) = colnames(xx)
  yy1 = yy1[, c(5,10:12)]
  
  yy2 = res[, c(1:10)]
  yy2 = yy2[, c(4:6,10)]
  ggs_intersect = intersect(rownames(yy1), rownames(yy2))
  
  yy1 = yy1[match(ggs_intersect, rownames(yy1)), ]
  yy2 = yy2[match(ggs_intersect, rownames(yy2)), ]
  
  yy1 = apply(yy1, 1, median)
  yy3 = apply(yy2, 1, median)
  
  ## test quantile normalization
  yy3 = preprocessCore::normalize.quantiles.use.target(as.matrix(yy2), target = as.vector(yy1))
  colnames(yy3) = colnames(yy2)
  rownames(yy3) = rownames(yy2)
  res = yy3
  
  yy3 = apply(yy3, 1, median)
  
  plot(yy1, yy3, cex = 0.5)
  abline(0, 1, lwd = 1.5, col = 'red')
  
  fit = lm(yy1 ~ yy3)
  
  yy4 = as.matrix(res[, c(1:10)]) * coefficients(fit)[2] + coefficients(fit)[1]
  
}

#saveRDS(genes2keep, file = paste0(RdataDir, 'genes_to_keep_RNAseq_vs_RPEChoroid.microarray.rds'))
yy = data.frame(xx0[match(genes2keep, rownames(xx0)), ], 
                res[match(genes2keep, rownames(res)), c(1:10)])

#yy = yy[, grep('GSM738', colnames(yy), invert = TRUE)]

pca2save <- prcomp(t(as.matrix(yy)), scale = FALSE)
#pca2save = pca2save[, -c(1:2)]
#pca2save = data.frame(pca2save$x, design0)
pca2save = data.frame(pca2save$x)
pca2save$name = rownames(pca2save)
#pca2save$age = factor(pca2save$age)

library(ggplot2)
ggplot(data=pca2save, aes(PC1, PC2, label = name))  + 
  geom_point(size=3) + 
  geom_text(hjust = 0.3, nudge_y = 0.3, size=2.5)

ggsave(paste0(resDir, '/RPE_RNAseq_vs_MA_PCA_age.genes.188.pdf'), 
       width=8, height = 4)



##########################################
# Combat
##########################################
Test_intergration_combat = FALSE

if(Test_intergration_combat){
  
  tmm = as.matrix(yy)
  
  Test_QN = FALSE
  if(Test_QN){
    yy2 = tmm[, grep('GSM', colnames(tmm))]
    yy0 = apply(tmm[, grep('syn_', colnames(tmm))], 1, mean)
    yy3 = preprocessCore::normalize.quantiles.use.target(as.matrix(yy2), target = as.vector(yy0))
    colnames(yy3) = colnames(yy2)
    rownames(yy3) = rownames(yy2)
    
    tmm = data.frame(tmm[, c(1:23)], yy3)
    tmm = as.matrix(tmm)
    
    
  }
  
  bc = as.factor(c(rep('RNAseq', nrow(design0)), rep('microarray', nrow(design))))
  conds = c('30', '40', '50', '60', '70', '80', '90', '60', '60', '70', '70', '70', '90',
            '0.3', 's10','s20', 's30','s40', 's50', 's60', 's70', 's80', 's90',
            '0.3', '0.3', '0.3', '70', '70', '70', '0.3', '0.3', '0.3', '70')
  mod = model.matrix(~ as.factor(conds))
  
  # if specify ref.batch, the parameters will be estimated from the ref, inapprioate here, 
  # because there is no better batche other others 
  #ref.batch = '2021S'# 2021S as reference is better for some reasons (NOT USED here)    
  fpm.bc = ComBat(dat=as.matrix(tmm), batch=bc, mod=mod, par.prior=FALSE, ref.batch = NULL) 
  
  pca2save <- prcomp(t(as.matrix(fpm.bc)), scale = FALSE)
  #pca2save = pca2save[, -c(1:2)]
  #pca2save = data.frame(pca2save$x, design0)
  pca2save = data.frame(pca2save$x)
  pca2save$name = rownames(pca2save)
  pca2save$age = factor(conds)
  pca2save$batch = bc
  
  library(ggplot2)
  ggplot(data=pca2save, aes(PC1, PC2, label = name, color = age, shape = batch))  + 
    geom_point(size=3) + 
    geom_text(hjust = 0.3, nudge_y = 0.3, size=2.5)
  
  ggsave(paste0(resDir, "/matureSamples_batchCorrect_before_",  version.analysis, ".pdf"), 
         width = 16, height = 14)
  
    
}


########################################################
########################################################
# Section V: Test quadratic programming to assign age
# 
########################################################
########################################################

##########################################
# import first the cellline RNA-seq data and process from Birgit
##########################################
require(DESeq2)

counts = readRDS('/groups/tanaka/Collaborations/Jingkui-Birgit/youngD0.rds')
counts2 = readRDS('/groups/tanaka/Collaborations/Jingkui-Birgit/oldD0.rds')
counts = data.frame(counts, counts2[match(rownames(counts), rownames(counts2)), ])

design = data.frame(age = rep(c('young', 'old'),  each = 4))
dds <- DESeqDataSetFromMatrix(as.matrix(counts), DataFrame(design), design = ~ age)

ss = rowSums(counts(dds))

hist(log10(ss), breaks = 100);abline(v = log10(20), lwd = 2.0, col = 'red')
cat(length(which(ss>20)), ' gene selected \n')
dds = dds[which(ss>20), ]

dds = estimateSizeFactors(dds)

saveRDS(dds, file = paste0(RdataDir, 'RNAseq_RPE.cellLine_Birgit.rds'))

##########################################
# test quadratic programming
##########################################
require(tibble)
require(tidyr)
require(dplyr)
library(gridExtra)
library(grid)
library(lattice)
library(ggpubr)
library(preprocessCore)

load(file = paste0(RdataDir, 'dds_design_ageGenes_RPE_RNAseq_Bulter_2021.Rdata'))
load(file = paste0(RdataDir,  'Bulter_2021_RPE_ageGenes_pval0.01_prediction.allTimepoints_design.Rdata'))

nb_genes_toUse = 50

Use_selecteGenes = TRUE
if(Use_selecteGenes){
  genes2keep = readRDS(file = paste0(RdataDir, 
                                     'genes_to_keep_RNAseq_vs_RPEChoroid.microarray_RPE.fetal.age.microarray.rds'))
  
  refs = xx0[match(genes2keep, rownames(xx0)), grep('syn_', colnames(xx0))]
  pvals = cpm$pval[match(rownames(refs), rownames(cpm))]
  refs$pvals = pvals
  refs = refs[order(refs$pvals), ]
  
  #refs = refs[c(1:50), c(1:10)]
  refs = refs[c(1:nb_genes_toUse), c(1:10)]
  
}

dds = readRDS(paste0(RdataDir, 'RNAseq_RPE.cellLine_Birgit.rds'))
xx = fpm(dds)
xx = as.matrix(log1p(xx))

Add_microarray = FALSE
if(Add_microarray){
  
  load(file = paste0(RdataDir, 'design_probAnnot_expr_gene.sample.matrix.Rdata'))
  yy1 = xx
    
  yy2 = res[, c(1:10)]
  #yy2 = yy2[, c(4:6,10)]
  ggs_intersect = intersect(rownames(yy1), rownames(yy2))
  
  yy1 = yy1[match(ggs_intersect, rownames(yy1)), ]
  yy2 = yy2[match(ggs_intersect, rownames(yy2)), ]
  
  xx = yy1
  yy1 = apply(yy1, 1, mean)
  yy2 = log1p(2^yy2)
  #yy3 = apply(yy2, 1, median)
  
  ## test quantile normalization
  yy3 = preprocessCore::normalize.quantiles.use.target(as.matrix(yy2), target = as.vector(yy1))
  colnames(yy3) = colnames(yy2)
  rownames(yy3) = rownames(yy2)
  
  res = yy3
  
  samples = cbind(apply(xx[, c(1:4)], 1, mean), apply(xx[, -c(1:4)], 1, mean))
  samples = cbind(samples, apply(res[, grep('w_', colnames(res))], 1, mean))
  samples = cbind(samples, apply(res[, grep('w_', colnames(res), invert = TRUE)], 1, mean))
  
  colnames(samples) = c('young', 'old', 'fetal.ma', 'old.ma')
  
}else{
  
  samples = xx
  #colnames(samples) = colnames(xx0)[1:13]
  #samples = samples[, order(design$age)]
  #samples = t(rbind(apply(xx[, c(1:4)], 1, mean), apply(xx[, -c(1:4)], 1, mean)))
  #colnames(samples) = c('young', 'old')
  
}

ggs_intersect = intersect(rownames(refs), rownames(samples))

samples = samples[match(ggs_intersect, rownames(samples)), ]
refs = refs[match(ggs_intersect, rownames(refs)), ]

source('functions_quadratic_programming.R')

samples = as.matrix(samples)
refs = as.matrix(log1p(2^refs - 2^-4))

given.cell.typs <- colnames(refs)

identity.matx <- data.frame()
for (i in 1:ncol(samples))
{
  # i = 1
  identity<-c()
  quad.rslt <- quad.prog.calc(i, refs, samples, force.eq = 0)
  
  QP <- quad.rslt[[1]]
  Error <- quad.rslt[[3]]
  
  identity <- c(colnames(samples)[i], QP$solution, QP$Lagrangian[1],Error)
  
  if (nrow(identity.matx) == 0) {
    identity.matx <- as.data.frame(t(as.matrix(identity)))
  } else {
    
    identity.matx <- rbind(identity.matx, as.data.frame(t(as.matrix(identity))))
  }
}

col.frx.names <- paste("frxn_", given.cell.typs, sep = "")
colnames(identity.matx)<-c("cell_name", col.frx.names,"Lagrangian","Error")

cat(length(ggs_intersect), ' overlapping genes used \n')
print(identity.matx$Error)

keep = t(identity.matx[, c(1:(1+length(given.cell.typs)))])
colnames(keep) = keep[1,]
keep = keep[-1, ]
keep = data.frame(keep)
keep2 = apply(keep, 2, function(x){x = as.numeric(x); x[which(x < 0.001)] = 0; x})
rownames(keep2) = rownames(keep)
keep = data.frame(keep2);rm(keep2)

keep$age = rownames(keep)
keep$age = gsub('frxn_syn_', '', keep$age)

age_levels = colnames(keep)

as_tibble(keep) %>%
  tidyr::gather(samples, probs, 1:(ncol(keep)-1)) %>%
  #ggplot(aes(x= factor(samples, levels = c('young', 'old')), y = probs, fill=age)) +
  ggplot(aes(x= factor(samples, levels = age_levels), y = probs, fill=age)) +
  geom_bar(stat="identity") +
  scale_fill_brewer(palette="Paired") + 
  labs(x = "", y = 'Probility') +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, angle = 90), 
        axis.text.y = element_text(size = 12))

ggsave(paste0(resDir, "/res_quadraticProgramming_cellLine_ageEstimate_geneNumber_", nrow(refs), ".pdf"), 
       width = 6, height = 4)
