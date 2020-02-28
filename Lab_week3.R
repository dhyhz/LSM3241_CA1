if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install()

lst = c("GEOquery", "affy", "limma", "huex10sttranscriptcluster.db", "org.Hs.eg.db", "oligo")

for(i in lst) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(i, version = "3.10")
}

library(GEOquery)
library(affy)
library(limma)
library(huex10sttranscriptcluster.db)
library(org.Hs.eg.db)
library(oligo)
install.packages(knitr)
library(knitr)
library(RColorBrewer)

gse <- getGEO('GSE143150', GSEMatrix=FALSE)

names(GSMList(gse))

filePaths <- getGEOSuppFiles('GSE143150')

file <- list.files('data/')
file
data <- read.celfiles(paste0('data/',file))
data

gsm <- GSMList(gse)[[1]]

names(Meta(gse))
names(Meta(gsm))

names(Meta(gsm))[!(names(Meta(gsm)) %in% names(Meta(gse)))]
Meta(gsm)[!(names(Meta(gsm)) %in% names(Meta(gse)))]

genotype <- function(gsm) { 
  Meta(gsm)[['characteristics_ch1']][2]
}
sapply(GSMList(gse),genotype)

treatment <- function(gsm){
  Meta(gsm)[['characteristics_ch1']][3]
}
sapply(GSMList(gse),treatment)

pd <- data.frame(genotype=as.factor(sapply(GSMList(gse),genotype)), treatment=as.factor(sapply(GSMList(gse),treatment)))
pd

pd$genotype <- as.factor(pd$genotype)
pd$treatment <- as.factor(pd$treatment)
levels(pd$genotype) <- c("KO","WT")
levels(pd$treatment) <- c("DM", "FM")
kable(pd)

files <- list.files("./data", pattern = ".CEL.gz")
celfiles <- paste0('data/',files)
celfiles
affydata <- read.celfiles(celfiles, phenoData = new("AnnotatedDataFrame",pd))
phenoData(affydata)

eset <- oligo::rma(affydata)
plotDensity(exprs(eset),xlab='log intensity',main="feature level densities after RMA",lwd=2)

pData(eset)

a <- paste(eset$genotype, eset$treatment)
a <- factor(a)
a
model.pre <- model.matrix( ~ 0 + a)
model.pre
## rename the model columns to correspond to the different growth conditions.
colnames(model.pre) <- c("DM_KO","FM_KO","DM_WT","FM_WT")

contrasts.pre <- makeContrasts((FM_WT - DM_WT),(FM_KO - DM_KO) , (FM_KO - FM_WT), (DM_KO - DM_WT), levels=model.pre)
contrasts.pre

fit.pre <- lmFit(eset, model.pre)
results = decideTests(fit.pre)
vennDiagram(results)
fitted.contrast.pre <- contrasts.fit(fit.pre,contrasts.pre)
fitted.ebayes.pre <- eBayes(fitted.contrast.pre)


## 1st Pairwise:  Full Media vs Diff Media on JMJD2B WT
topTable(fitted.ebayes.pre, coef = 1)

ps2.pre1 <- topTable(fitted.ebayes.pre,coef=1,number=Inf,p.value = 0.05,lfc=2)
ps2_up.pre1 <- rownames(ps2.pre1[ps2.pre1$logFC>2,])
AnnotationDbi::select(huex10sttranscriptcluster.db,ps2_up.pre1,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")
write.csv(AnnotationDbi::select(huex10sttranscriptcluster.db,ps2_up.pre1,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID"), "genespre1up.csv")

ps2.pre1 <- topTable(fitted.ebayes.pre,coef=1,number=Inf,p.value = 0.05,lfc=2)
ps2_up.pre1 <- rownames(ps2.pre1[ps2.pre1$logFC<-2,])
AnnotationDbi::select(huex10sttranscriptcluster.db,ps2_up.pre1,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")
write.csv(AnnotationDbi::select(huex10sttranscriptcluster.db,ps2_up.pre1,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID"), "genespre1down.csv")

ps2_up.pre1 <- rownames(ps2.pre1)
AnnotationDbi::select(huex10sttranscriptcluster.db,ps2_up.pre1,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")
write.csv(AnnotationDbi::select(huex10sttranscriptcluster.db,ps2_up.pre1,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID"), "genespre1.csv")

volcanoplot(fitted.ebayes.pre)

interesting_genes.pre1 <- topTable(fitted.ebayes.pre,coef = 1,number=Inf,p.value = 0.05,lfc=2)
volcanoplot(fitted.ebayes.pre, coef=1, main=sprintf("%d features pass our cutoffs",nrow(interesting_genes.pre1)))
points(interesting_genes.pre1[['logFC']],-log10(interesting_genes.pre1[['P.Value']]),col='red')

eset_of_interest.pre1 <- eset[rownames(interesting_genes.pre1),]
heatmap(exprs(eset_of_interest.pre1),
        labCol=eset$treatment,labRow=NA,
        col       = rev(brewer.pal(10, "RdBu")),
        distfun   = function(x) as.dist(1-cor(t(x))))

#############################################################


## 2nd Pairwise: Full Media vs Diff Media on JMJD2B KO 
topTable(fitted.ebayes.pre, coef = 2)

ps2.pre2 <- topTable(fitted.ebayes.pre,coef=2,number=Inf,p.value = 0.05,lfc=2)
ps2_up.pre2 <- rownames(ps2.pre2[ps2.pre2$logFC>2,])
AnnotationDbi::select(huex10sttranscriptcluster.db,ps2_up.pre2,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")
write.csv(AnnotationDbi::select(huex10sttranscriptcluster.db,ps2_up.pre2,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID"), "genespre2up.csv")

ps2.pre2 <- topTable(fitted.ebayes.pre,coef=2,number=Inf,p.value = 0.05,lfc=2)
ps2_up.pre2 <- rownames(ps2.pre2[ps2.pre2$logFC<-2,])
AnnotationDbi::select(huex10sttranscriptcluster.db,ps2_up.pre2,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")
write.csv(AnnotationDbi::select(huex10sttranscriptcluster.db,ps2_up.pre2,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID"), "genespre2down.csv")

ps2.pre2 <- topTable(fitted.ebayes.pre,coef=2,number=Inf,p.value = 0.05,lfc=2)
ps2_up.pre2 <- rownames(ps2.pre2)
AnnotationDbi::select(huex10sttranscriptcluster.db,ps2_up.pre2,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")
write.csv(AnnotationDbi::select(huex10sttranscriptcluster.db,ps2_up.pre2,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID"), "genespre2.csv")

volcanoplot(fitted.ebayes.pre)

interesting_genes.pre2 <- topTable(fitted.ebayes.pre, coef=2, number=Inf,p.value = 0.05,lfc=2)
volcanoplot(fitted.ebayes.pre, coef=2, main=sprintf("%d features pass our cutoffs",nrow(interesting_genes.pre2)))
points(interesting_genes.pre2[['logFC']],-log10(interesting_genes.pre2[['P.Value']]),col='red')

eset_of_interest.pre2 <- eset[rownames(interesting_genes.pre2),]
heatmap(exprs(eset_of_interest2.pre),
        labCol=eset$treatment,labRow=NA,
        col       = rev(brewer.pal(10, "RdBu")),
        distfun   = function(x) as.dist(1-cor(t(x))))

#############################################################


## 3RD Pairwise:  Full Media vs Diff Media on JMJD2B WT
topTable(fitted.ebayes.pre, coef = 3)

ps2.pre3 <- topTable(fitted.ebayes.pre,coef=3,number=Inf,p.value = 0.05,lfc=1)
ps2_up.pre3 <- rownames(ps2.pre3)
AnnotationDbi::select(huex10sttranscriptcluster.db,ps2_up.pre3,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")

volcanoplot(fitted.ebayes.pre)

interesting_genes.pre3 <- topTable(fitted.ebayes.pre,coef = 3,number=Inf,p.value = 0.05,lfc=2)
volcanoplot(fitted.ebayes.pre, coef=3, main=sprintf("%d features pass our cutoffs",nrow(interesting_genes.pre3)))
points(interesting_genes.pre3[['logFC']],-log10(interesting_genes.pre3[['P.Value']]),col='red')

eset_of_interest.pre3 <- eset[rownames(interesting_genes.pre3),]
heatmap(exprs(eset_of_interest3.pre),
        labCol=eset$genotype,labRow=NA,
        col       = rev(brewer.pal(10, "RdBu")),
        distfun   = function(x) as.dist(1-cor(t(x))))

#############################################################


## 4TH Pairwise:  Full Media vs Diff Media on JMJD2B WT
topTable(fitted.ebayes.pre, coef = 4)

ps2.pre4 <- topTable(fitted.ebayes.pre,coef=4,number=Inf,p.value = 0.05,lfc=1)
ps2_up.pre4 <- rownames(ps2.pre4)
AnnotationDbi::select(huex10sttranscriptcluster.db,ps2_up.pre4,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")

volcanoplot(fitted.ebayes.pre)

interesting_genes.pre4 <- topTable(fitted.ebayes.pre4,coef = 4,number=Inf,p.value = 0.05,lfc=2)
volcanoplot(fitted.ebayes.pre, coef=4, main=sprintf("%d features pass our cutoffs",nrow(interesting_genes.pre4)))
points(interesting_genes.pre4[['logFC']],-log10(interesting_genes.pre4[['P.Value']]),col='red')

eset_of_interest.pre4 <- eset[rownames(interesting_genes.pre4),]
heatmap(exprs(eset_of_interest.pre4),
        labCol=eset$genotype,labRow=NA,
        col       = rev(brewer.pal(10, "RdBu")),
        distfun   = function(x) as.dist(1-cor(t(x))))

#############################################################



## Final model
model.final <- model.matrix( ~ 0 + a)
model.final
colnames(model.final)<-c("DM_KO","FM_KO","DM_WT","FM_WT")

contrasts.final <- makeContrasts((FM_KO + FM_WT) - (DM_KO + DM_WT), levels=model.final)
contrasts.final

fit.final <- lmFit(eset, model.final)
fitted.contrast.final <- contrasts.fit(fit.final,contrasts.final)

fitted.ebayes.final <- eBayes(fitted.contrast.final)
topTable(fitted.ebayes.final)

ps2.final <- topTable(fitted.ebayes.final,number=Inf,p.value = 0.05,lfc=2)
ps2_up.final <- rownames(ps2.final[ps2.final$logFC<-2,])
write.csv(AnnotationDbi::select(huex10sttranscriptcluster.db,ps2_up.final,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID"), "genesfinaldown.csv")
AnnotationDbi::select(huex10sttranscriptcluster.db,ps2_up.final,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")

ps2.final <- topTable(fitted.ebayes.final,number=Inf,p.value = 0.05,lfc=2)
ps2_up.final <- rownames(ps2.final[ps2.final$logFC>2,])
write.csv(AnnotationDbi::select(huex10sttranscriptcluster.db,ps2_up.final,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID"), "genesfinalup.csv")
AnnotationDbi::select(huex10sttranscriptcluster.db,ps2_up.final,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")

ps2.final <- topTable(fitted.ebayes.final,number=Inf,p.value = 0.05,lfc=2)
ps2_up.final <- rownames(ps2.final)
write.csv(AnnotationDbi::select(huex10sttranscriptcluster.db,ps2_up.final,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID"), "genesfinal.csv")
AnnotationDbi::select(huex10sttranscriptcluster.db,ps2_up.final,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")

volcanoplot(fitted.ebayes.final)

interesting_genes.final <- topTable(fitted.ebayes.final,number=Inf,p.value = 0.05,lfc=2)
volcanoplot(fitted.ebayes.final, main=sprintf("%d features pass our cutoffs",nrow(interesting_genes.final)))
points(interesting_genes.final[['logFC']],-log10(interesting_genes.final[['P.Value']]),col='red')

eset_of_interest.final <- eset[rownames(interesting_genes.final),]
heatmap(exprs(eset_of_interest.final),
        labCol=eset$treatment,labRow=NA,
        col       = rev(brewer.pal(10, "RdBu")),
        distfun   = function(x) as.dist(1-cor(t(x))))
