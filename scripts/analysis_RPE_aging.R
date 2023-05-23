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
# Section I : process micorarray data from GSE18811
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
#bg = mat[, grep('_BG', colnames(mat))]
#mat = mat[, grep('_BG', colnames(mat), invert = TRUE)]

#ex <- log2(ex)  # log2 transformation

# # box-and-whisker plot
# dev.new(width=3+ncol(gset)/6, height=5)
# par(mar=c(7,4,2,1))
# title <- paste ("GSE18811", "/", annotation(gset), sep ="")
# boxplot(ex, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)
# dev.off()
# 
# # expression value distribution plot
# par(mar=c(4,4,2,1))
# title <- paste ("GSE18811", "/", annotation(gset), " value distribution", sep ="")
# plotDensities(ex, main=title, legend=F)
# 
# # mean-variance trend
# ex <- na.omit(ex) # eliminate rows with NAs
# plotSA(lmFit(ex), main="Mean variance trend, GSE18811")
# 
# # UMAP plot (multi-dimensional scaling)
# ex <- ex[!duplicated(ex), ]  # remove duplicates
# ump <- umap(t(ex), n_neighbors = 13, random_state = 123)
# plot(ump$layout, main="UMAP plot, nbrs=13", xlab="", ylab="", pch=20, cex=1.5)
# library("maptools")  # point labels without overlaps
# pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)

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
# Section 1.5:
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

pdf(paste0(resDir, "/RNAseq_ageGenes_pval0.01_vs_Newman.2012.pdf"), height = 8, width =10)

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


########################################################
########################################################
# Section III: process RNA-seq data from Bulter_2021_RNAseq_GSE159435
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
require(ggplot2)

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
colnames(xx) = c('pval', 'beta', 'intercept')

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
# Section II: 
# test to integrate RNA-seq data and microarray 
# 
########################################################
########################################################

##########################################
# test Quantile normalization 
##########################################
Test_intergration_QN = FALSE
if(Test_intergration_QN){
  xx1 = preprocessCore::normalize.quantiles.use.target(as.matrix(xx2), target = as.vector(xx0[,1]))
  colnames(xx1) = colnames(xx2)
  rownames(xx1) = rownames(xx2)
  
  tmm = as.matrix(cbind(xx0, xx1))
  #tmm.qn = preprocessCore::normalize.quantiles(tmm)
  #colnames(tmm.qn) = colnames(tmm)
  #rownames(tmm.qn) = rownames(tmm)
  #tmm = tmm.qn
  #make.pca.plots(tmm, ntop = 3000)
  bc = as.factor(c(rep('RNAseq', nrow(design0)), rep('microarray', nrow(design))))
  conds = c('31', '40', '51', '67', '70.80', '86', '93', '61', '62', '70.80', '70.80', '70.80', '92',
            '0.3', '0.3', '0.3', '70.80', '70.80', '70.80', '0.3', '0.3', '0.3', '70.80')
  mod = model.matrix(~ as.factor(conds))
  
  # if specify ref.batch, the parameters will be estimated from the ref, inapprioate here, 
  # because there is no better batche other others 
  #ref.batch = '2021S'# 2021S as reference is better for some reasons (NOT USED here)    
  fpm.bc = ComBat(dat=as.matrix(tmm), batch=bc, mod=mod, par.prior=TRUE, ref.batch = NULL) 
  
  make.pca.plots(fpm.bc, pca.dim = c(1, 3), ntop = 500)
  #make.pca.plots(fpm.bc, ntop = 3000)
  ggsave(paste0(resDir, "/first_test_batchCorrect_RNAseq_microarray.pdf"), width = 10, height = 8)
  
}

load(file = paste0(RdataDir, 'dds_design_ageGenes_RPE_RNAseq_Bulter_2021.Rdata'))
design0 = design

load(file = paste0(RdataDir, 'design_probAnnot_expr_gene.sample.matrix.Rdata'))
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

xx0 = assay(vsd)
xx2 = res[, c(1:10)]

ggs = intersect(rownames(xx0), rownames(xx2))

xx0 = xx0[match(ggs, rownames(xx0)),]
xx2 = xx2[match(ggs, rownames(xx2)), ]

colnames(xx0) = paste0(design0$donor, '_', design0$age, '_',  design0$sex)


##########################################
# Combat
##########################################
Test_intergration_combat = FALSE

if(Test_intergration_combat){
  tmm = as.matrix(cbind(xx0, xx2))
  
  bc = as.factor(c(rep('RNAseq', nrow(design0)), rep('microarray', nrow(design))))
  conds = c('31', '40', '51', '67', '70.80', '86', '93', '61', '62', '70.80', '70.80', '70.80', '92',
            '0.3', '0.3', '0.3', '70.80', '70.80', '70.80', '0.3', '0.3', '0.3', '70.80')
  mod = model.matrix(~ as.factor(conds))
  
  # if specify ref.batch, the parameters will be estimated from the ref, inapprioate here, 
  # because there is no better batche other others 
  #ref.batch = '2021S'# 2021S as reference is better for some reasons (NOT USED here)    
  fpm.bc = ComBat(dat=as.matrix(tmm), batch=bc, mod=mod, par.prior=TRUE, ref.batch = NULL) 
  
  make.pca.plots(tmm, ntop = 3000)
  ggsave(paste0(resDir, "/matureSamples_batchCorrect_before_",  version.analysis, ".pdf"), 
         width = 16, height = 14)
  
  make.pca.plots(fpm.bc, ntop = 3000)
  ggsave(paste0(resDir, "/first_test_batchCorrect_RNAseq_microarray.pdf"), width = 16, height = 14)
  
}



