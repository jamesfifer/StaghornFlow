 source("http://bioconductor.org/biocLite.R")
 biocLite("DESeq2")

###conduct array quality metrics to detect and remove outliers

library(DESeq) 
library(affycoretools)
library(arrayQualityMetrics)
library(genefilter)
library(DESeq2)
library(Biobase)
library(pheatmap)
library(vegan)
library(rgl)
library(ape)
library(VennDiagram)
 BiocManager::install("apeglm")
library(apeglm)

#read in raw counts from rsem and do correction by length via tximport
 
 if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
 
 BiocManager::install("tximport")
 library(tximport)
 setwd('C:/Users/james/Documents/Documents/JAMES Broken Computer/Guam/Raymundo/Reanalysis/Ex_Situ')
# dir <- 'C:/Users/james/Documents/BOSTON/Ecological & evolutionary genomics/Final_Project'
 dir <- 'C:/Users/james/Documents/Documents/JAMES Broken Computer/Guam/Raymundo/Reanalysis/Ex_Situ'
 
  list.files(dir)
 #samples.txt should at minimum have the directory or file name headers (e.g. HFH_rep1)
 samples <- read.table(file.path(dir, "samples.txt"), header = FALSE)
 
 #Use this code if each file is still within its own folder
 #files <- file.path(dir, "rsem", samples$V2, paste0(samples$V2, ".RSEM.genes.results"))
 
 #make sure the column you call (V2) matches the first part of the file names in your directory.
 files <- file.path(dir, paste0(samples$V2, ".RSEM.genes.results"))
 files
 names(files) <- paste0(samples$V2)
 all(file.exists(files))
 txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
 head(txi.rsem$counts)
 txi.rsem
 #Note: there are two suggested ways of importing estimates for use with differential gene expression (DGE) methods.
 #The first method, which we show below for edgeR and for DESeq2, is to use the gene-level estimated counts from the
 #quantification tools, and additionally to use the transcript-level abundance estimates to calculate a gene-level offset
 #that corrects for changes to the average transcript length across samples. The code examples below accomplish these steps for you,
 #keeping track of appropriate matrices and calculating these offsets. For edgeR you need to assign a matrix to y$offset, but the function
 #DESeqDataSetFromTximport takes care of creation of the offset for you. Let's call this method "original counts and offset".
 
 #replace HFH with HFH, HFNH with HFNH, LFH with LFH, LFH with LFNH. 
 treat=c( "HFH", "HFH", "HFH" , "HFH", "HFNH","HFNH", "HFNH", "HFNH", "LFH","LFH","LFH","LFH","LFNH","LFNH","LFNH")
 colony=as.factor (c( "Col_b", "Col_a", "Col_c", "Col_d", "Col_a","Col_b", "Col_c", "Col_d", "Col_b","Col_a","Col_c","Col_d","Col_b","Col_a","Col_c"))
 g=data.frame(treat, colony)
 g
 colData<- g
 
 #for heat comp
 heat= treat=c( "Heat", "Heat", "Heat" , "Heat", "Ambient","Ambient", "Ambient", "Ambient", "Heat","Heat","Heat","Heat","Ambient","Ambient","Ambient")
flow=c( "HF", "HF", "HF" , "HF", "HF","HF", "HF", "HF", "LF","LF","LF","LF","LF","LF","LF")
colony=as.factor (c( "Col_b", "Col_a", "Col_c", "Col_d", "Col_a","Col_b", "Col_c", "Col_d", "Col_b","Col_a","Col_c","Col_d","Col_b","Col_a","Col_c"))
 g=data.frame(heat,flow, colony)
 g
 colData<- g
 
 
 head(txi.rsem$counts)
 length(txi.rsem$counts[,1])
 #27807 host
 #22495 sym
 
 #would need to define m as an integer, not sure how kosher it is to just integer these expected counts(txi.rsem$counts)
 m<-txi.rsem$counts
 #
 library(dplyr)
 library(tibble)  # for `rownames_to_column` and `column_to_rownames`
 

 
 storage.mode(m) = "integer"
 conditions=data.frame(treat)
 
 real=newCountDataSet(m,conditions) 
 real=estimateSizeFactors(real)
 plot(sort(sizeFactors(real))) 
 
 cds=estimateDispersions(real,method="blind")
 vsdBlind=varianceStabilizingTransformation(cds)
 
 v=setwd("arrayqmetrics/")
 
 arrayQualityMetrics(vsdBlind,intgroup=c("treat"), force=TRUE)
 # # # #look for outliers
 
 #####you only ever need to run the above code once. Outliers are decided at the beginning. 
 ## I like to close R and restart with packages etcw
 ###HFNH_rep_1 had only one star so is not considered as an outlier 
 
#load packages
 detach("package:DESeq", unload=TRUE)
 library(ggplot2)
 library(Biobase)
 library(pheatmap)
 library(vegan)
 library(rgl)
 library(ape)
 library(VennDiagram)
 library(dplyr)
 library(DESeq2) 
 library(pheatmap)
 library(ggrepel)
 library(tidyverse)



 #treat+genotype is what Michael Love uses, similar to  treat+batch to "control for batch effects"
 #This model was for comparisons between all four treatments (HFNH, HFH, LFNH, LFH)
 dds<-DESeqDataSetFromTximport(txi.rsem, colData, design=~treat+colony)
 
 
 #Thisnext  model was for looking for signature of heat stress
 
 #dds<-DESeqDataSetFromTximport(txi.rsem, colData, design=~heat+flow+colony+heat:flow)
 #adding the interaction effect does not actual control for the interaction so the above model was not used
 #adding the interaction effect and then looking at comparisons within heat for example
 #will change the main effect of heat so you are only looking at the first reference level of flow (LF in this case)
 #see vignette for more info http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
 
 dds<-DESeqDataSetFromTximport(txi.rsem, colData, design=~heat+flow+colony)
 

#deseq2 was not filtering out some low guys
keep <- rowSums(counts(dds) >= 10) >= 3
dds<-dds[keep,]



totalCounts=as.vector(colSums(m))
totalCounts

dev.off()
barplot(totalCounts, col=c("coral", "coral", "coral", "coral", "red", "red", "red", "red", "blue", "blue", "blue", "blue", "green", "green" ,"green"), ylab="raw counts")

totalCounts
#HFH  HFH  HFH  HFH  HFNH  HFNH  HFNH  HFNH  LFH  LFH  LFH  LFH  LFNH  LFNH  LFNH 

#2388512  937645  947137 1660043 3707369 1460498 1665146  222477  431148 1385240 2088889 1237827  336166  434172  220117


min(totalCounts) #220117
max(totalCounts)  # 3707369


#run deseq
dds<-DESeq(dds)


head(dds)
res<- results(dds)

#Look at dispersions plot
plotDispEsts(dds, main="Dispersion plot ")
colData$colony
colData$heat

################Pooled heat and flow comparisons #################
#using model dds<-DESeqDataSetFromTximport(txi.rsem, colData, design=~heat+flow+colony)
#heat comparison
colData$HvsNH<-factor(colData$heat, levels=c("Heat","Ambient"))
resHvsNH <- results(dds, contrast=c("heat","Heat","Ambient"))

table(resHvsNH$padj<0.1)
table(resHvsNH$padj<0.05)
table(resHvsNH$padj<0.01)

nrow(resHvsNH[resHvsNH$padj<0.05 & resHvsNH$log2FoldChange > 0 & !is.na(resHvsNH$padj),])
nrow(resHvsNH[resHvsNH$padj<0.05 & resHvsNH$log2FoldChange < 0 & !is.na(resHvsNH$padj),])
#518
#414

write.table(resHvsNH, file="HvsNH.txt", quote=F, sep="\t")

cd <- read.table("HvsNH.txt")

#Add annotations
ID2gene<-read.csv(file = "ID2gene.tab", sep="\t", row.names=2)
ID2gene$row.names<-NULL
newcd<-transform(merge(cd,ID2gene,by=0), row.names=Row.names, Row.names=NULL)
#annotated
nrow(newcd[newcd$padj<0.05 & newcd$log2FoldChange > 0 & !is.na(newcd$padj),])
nrow(newcd[newcd$padj<0.05 & newcd$log2FoldChange < 0 & !is.na(newcd$padj),])
#281
#277


write.table(newcd, file="HvsNH.txt")


MAttPlot <- function(df) {
  df$dotcol <- ifelse(df$log2FoldChange > 0 & df$padj < 0.1, 'darkorange',
                      ifelse(df$log2FoldChange < 0 & df$padj < 0.1, 'cyan3', 'black'))
  df$baseMean <- log2(df$baseMean)
  print(head(df))
  gg <- ggplot(df, aes(baseMean, log2FoldChange)) +
    geom_point(size = .3, color = df$dotcol) +
    theme_bw() +
    theme(panel.grid = element_blank())
  print(gg)
}

MAttPlot(cd)

##########make the GO table for MWU
head(cd)
cd$isogroup=row.names(cd)
library(dplyr)
go_input_HvsNH = cd %>%
  mutate(mutated_p = -log(pvalue)) %>%
  mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
  select(isogroup, mutated_p_updown) %>%
  na.omit()

nrow(go_input_HvsNH)
head(go_input_HvsNH)
nrow(go_input_HvsNH)
colnames(go_input_HvsNH) <- c("gene", "pval")
head(go_input_HvsNH)
write.csv(go_input_HvsNH, file="HvsNH_GO.csv", quote=F, row.names=FALSE)


colData$flow
##############       FLOW        ################################### FLOW #########
# using model dds<-DESeqDataSetFromTximport(txi.rsem, colData, design=~heat+flow+colony)
colData$HFvsLF<-factor(colData$flow, levels=c("HF","LF"))
resHFvsLF <- results(dds, contrast=c("flow","HF","LF"))

table(resHFvsLF$padj<0.1)
table(resHFvsLF$padj<0.05)
table(resHFvsLF$padj<0.01)

nrow(resHFvsLF[resHFvsLF$padj<0.05 & resHFvsLF$log2FoldChange > 0 & !is.na(resHFvsLF$padj),])
nrow(resHFvsLF[resHFvsLF$padj<0.05 & resHFvsLF$log2FoldChange < 0 & !is.na(resHFvsLF$padj),])
#518
#414

write.table(resHFvsLF, file="HFvsLF.txt", quote=F, sep="\t")

cd <- read.table("HFvsLF.txt")

#Add annotations
ID2gene<-read.csv(file = "ID2gene.tab", sep="\t", row.names=2)
ID2gene$row.names<-NULL
newcd<-transform(merge(cd,ID2gene,by=0), row.names=Row.names, Row.names=NULL)
#annotated
nrow(newcd[newcd$padj<0.05 & newcd$log2FoldChange > 0 & !is.na(newcd$padj),])
nrow(newcd[newcd$padj<0.05 & newcd$log2FoldChange < 0 & !is.na(newcd$padj),])
#281
#277

write.table(newcd, file="HFvsLF.txt")


MAttPlot <- function(df) {
  df$dotcol <- ifelse(df$log2FoldChange > 0 & df$padj < 0.1, 'darkorange',
                      ifelse(df$log2FoldChange < 0 & df$padj < 0.1, 'cyan3', 'black'))
  df$baseMean <- log2(df$baseMean)
  print(head(df))
  gg <- ggplot(df, aes(baseMean, log2FoldChange)) +
    geom_point(size = .3, color = df$dotcol) +
    theme_bw() +
    theme(panel.grid = element_blank())
  print(gg)
}

MAttPlot(cd)

##########make the GO table for MWU
head(cd)
cd$isogroup=row.names(cd)
library(dplyr)
go_input_HFvsLF = cd %>%
  mutate(mutated_p = -log(pvalue)) %>%
  mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
  select(isogroup, mutated_p_updown) %>%
  na.omit()

nrow(go_input_HFvsLF)
head(go_input_HFvsLF)
nrow(go_input_HFvsLF)
colnames(go_input_HFvsLF) <- c("gene", "pval")
head(go_input_HFvsLF)
write.csv(go_input_HFvsLF, file="HFvsLF_GO.csv", quote=F, row.names=FALSE)

###############interaction##############################
#interaction comp (from model flow+heat+colony+flow:heat)
#this is not the best way to look for interactive effect between temp and flow, we are limiting our
#comparisons here, so this was not actually used in the analysis 

res=results(dds, name="heatHeat.flowLF")
res=results(dds, name="flowLF.heatAmbient")

#measures is the flow effect different across temperatures?
#log2(LF/HF(Heat)) - log2(LF/HF(Ambient))
#when the difference between lf and hf is greater under heat than under ambient
#it is positive 


res
table(res$padj<0.1)
#FALSE  TRUE 
#21723     5
table(res$padj<0.05)
# TRUE 
#   5
table(res$padj<0.01)
#4

write.table(res, file="heatHeat.flowLF", quote=F, sep="\t")
cd <- read.table("heatHeat.flowLF")
#Add annotations
ID2gene<-read.csv(file = "ID2gene.tab", sep="\t", row.names=2)
ID2gene$row.names<-NULL
#newcd<-merge(cd, ID2gene)

newcd<-transform(merge(cd,ID2gene,by=0, row.names=Row.names, Row.names=NULL))

write.table(newcd, file="heatHeat.flowLF.txt")

############################# A Comparisons  #######################################################
#################### A vs B 
colData$AvsB<-factor(colData$treat, levels=c("HFH","HFNH"))
##second term is the "control"
resAvsB <- results(dds, contrast=c("treat","HFH","HFNH"))
#how many FDR < 10%?

table(resAvsB$padj<0.1)
table(resAvsB$padj<0.05)
table(resAvsB$padj<0.01)
# 0.1=0
# 0.05=0
# 0.01=0
summary(resAvsB)
colnames(resAvsB)
nrow(resAvsB[resAvsB$padj<0.05 & !is.na(resAvsB$padj),])  # Num significantly differentially expressed genes excluding the no/low count genes   #228

plotMA(resAvsB, main="HFH vs HFNH")
plotMA(resAvsB, main="HFH vs HFNH + y-lim", ylim=c(-2,2))



results <- as.data.frame(resAvsB)
head(results)

nrow(resAvsB[resAvsB$padj<0.05 & resAvsB$log2FoldChange > 0 & !is.na(resAvsB$padj),])
nrow(resAvsB[resAvsB$padj<0.05 & resAvsB$log2FoldChange < 0 & !is.na(resAvsB$padj),])
#165
#181
test=data.frame(resAvsB[resAvsB$padj<0.05 & resAvsB$log2FoldChange > 0 & !is.na(resAvsB$padj),])
head(test)


write.table(resAvsB, file="AvsB.txt", quote=F, sep="\t")

cd <- read.table("AvsB.txt")

#Add annotations
ID2gene<-read.csv(file = "ID2gene.tab", sep="\t", row.names=2)
ID2gene$row.names<-NULL
newcd<-transform(merge(cd,ID2gene,by=0), row.names=Row.names, Row.names=NULL)
#annotated
nrow(newcd[newcd$padj<0.05 & newcd$log2FoldChange > 0 & !is.na(newcd$padj),])
nrow(newcd[newcd$padj<0.05 & newcd$log2FoldChange < 0 & !is.na(newcd$padj),])
#76
#145

write.table(newcd, file="AvsB.txt")


MAttPlot <- function(df) {
	df$dotcol <- ifelse(df$log2FoldChange > 0 & df$padj < 0.1, 'darkorange',
	ifelse(df$log2FoldChange < 0 & df$padj < 0.1, 'cyan3', 'black'))
	df$baseMean <- log2(df$baseMean)
	print(head(df))
	gg <- ggplot(df, aes(baseMean, log2FoldChange)) +
	geom_point(size = .3, color = df$dotcol) +
	theme_bw() +
	theme(panel.grid = element_blank())
	print(gg)
}

MAttPlot(cd)

##########make the GO table for MWU
head(cd)
cd$isogroup=row.names(cd)
library(dplyr)
go_input_AvsB = cd %>%
  mutate(mutated_p = -log(pvalue)) %>%
  mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
  select(isogroup, mutated_p_updown) %>%
  na.omit()

nrow(go_input_AvsB)
head(go_input_AvsB)
nrow(go_input_AvsB)
colnames(go_input_AvsB) <- c("gene", "pval")
head(go_input_AvsB)
write.csv(go_input_AvsB, file="AvsB_GO.csv", quote=F, row.names=FALSE)

#################### A vs C
colData$AvsC<-factor(colData$treat, levels=c("HFH","LFH"))
##second term is the "control"
resAvsC <- results(dds, contrast=c("treat","HFH","LFH"))
#how many FDR < 10%?
table(resAvsC$padj<0.1)
table(resAvsC$padj<0.05)
table(resAvsC$padj<0.01)
# 0.1= 1
# 0.05=0
# 0.01=0
summary(resAvsC)

nrow(resAvsC[resAvsC$padj<0.05 & !is.na(resAvsC$padj),])  # Num significantly differentially expressed genes excluding the no/low count genes   #228

plotMA(resAvsC, main="HFH vs LFH")
plotMA(resAvsC, main="HFH vs LFH + y-lim", ylim=c(-2,2))

results <- as.data.frame(resAvsC)
head(results)

nrow(resAvsC[resAvsC$padj<0.05 & resAvsC$log2FoldChange > 0 & !is.na(resAvsC$padj),])
nrow(resAvsC[resAvsC$padj<0.05 & resAvsC$log2FoldChange < 0 & !is.na(resAvsC$padj),])
#UP in A=0
#DOWN in A=3

write.table(resAvsC, file="AvsC.txt", quote=F, sep="\t")
#Add annotations
cd <- read.table("AvsC.txt")
ID2gene<-read.csv(file = "ID2gene.tab", sep="\t", row.names=2)
ID2gene$row.names<-NULL
newcd<-transform(merge(cd,ID2gene,by=0), row.names=Row.names, Row.names=NULL)
nrow(newcd[newcd$padj<0.05 & newcd$log2FoldChange > 0 & !is.na(newcd$padj),])
nrow(newcd[newcd$padj<0.05 & newcd$log2FoldChange < 0 & !is.na(newcd$padj),])
#0
#2
write.table(newcd, file="AvsC.txt")


MAttPlot <- function(df) {
  df$dotcol <- ifelse(df$log2FoldChange > 0 & df$padj < 0.1, 'darkorange',
                      ifelse(df$log2FoldChange < 0 & df$padj < 0.1, 'cyan3', 'black'))
  df$baseMean <- log2(df$baseMean)
  print(head(df))
  gg <- ggplot(df, aes(baseMean, log2FoldChange)) +
    geom_point(size = .3, color = df$dotcol) +
    theme_bw() +
    theme(panel.grid = element_blank())
  print(gg)
}

MAttPlot(cd)

##########make the GO table for MWU
head(cd)
cd$isogroup=row.names(cd)
library(dplyr)
go_input_AvsC = cd %>%
  mutate(mutated_p = -log(pvalue)) %>%
  mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
  select(isogroup, mutated_p_updown) %>%
  na.omit()

nrow(go_input_AvsC)
head(go_input_AvsC)
nrow(go_input_AvsC)
colnames(go_input_AvsC) <- c("gene", "pval")
head(go_input_AvsC)
write.csv(go_input_AvsC, file="AvsC_GO.csv", quote=F, row.names=FALSE)


#################### B vs D  #The better flow comparison# (no heated samples)
colData$BvsD<-factor(colData$treat, levels=c("HFNH","LFNH"))
##second term is the "control"
#resBvsD <- results(dds, contrast=c("treat","HFNH","LFNH"))
resLFC <- lfcShrink(dds, contrast=c("treat","HFNH","LFNH"), type="normal")
#how many FDR < 10%?
table(resLFC$padj<0.1)
table(resLFC$padj<0.05)
table(resLFC$padj<0.01)
# 0.1= 2851
# 0.05=1513
# 0.01=294
summary(resBvsD)


nrow(resBvsD[resBvsD$padj<0.05 & !is.na(resBvsD$padj),])  # Num significantly differentially expressed genes excluding the no/low count genes   
#3576

plotMA(resLFC, main="HFNH vs LFNH")
plotMA(resBvsD, main="HFNH vs LFNH + y-lim", ylim=c(-2,2))

cd=as.data.frame(resLFC)

nrow(resLFC[resLFC$padj<0.05 & resLFC$log2FoldChange > 0 & !is.na(resLFC$padj),])
nrow(resLFC[resLFC$padj<0.05 & resLFC$log2FoldChange < 0 & !is.na(resLFC$padj),])
#UP in B 3460
#DOWN in B 116

write.table(resLFC, file="BvsD.txt", quote=F, sep="\t")

cd <- read.table("BvsD.txt")
#Add annotations
ID2gene<-read.csv(file = "ID2gene.tab", sep="\t", row.names=2)
ID2gene$row.names<-NULL
newcd<-transform(merge(cd,ID2gene,by=0), row.names=Row.names, Row.names=NULL)
nrow(newcd[newcd$padj<0.05 & newcd$log2FoldChange > 0 & !is.na(newcd$padj),])
nrow(newcd[newcd$padj<0.05 & newcd$log2FoldChange < 0 & !is.na(newcd$padj),])
#3129
#42

write.table(newcd, file="BvsD.txt")

head(cd)

MAttPlot <- function(df) {
  df$dotcol <- ifelse(df$log2FoldChange > 0 & df$padj < 0.1, 'darkorange',
                      ifelse(df$log2FoldChange < 0 & df$padj < 0.1, 'cyan3', 'black'))
  df$baseMean <- log2(df$baseMean)
  print(head(df))
  gg <- ggplot(df, aes(baseMean, log2FoldChange)) +
    geom_point(size = .3, color = df$dotcol) +
    theme_bw() +
    theme(panel.grid = element_blank())
  print(gg)
}

MAttPlot(cd)

##########make the GO table for MWU
head(cd)
cd$isogroup=row.names(cd)
library(dplyr)
go_input_BvsD = cd %>%
  mutate(mutated_p = -log(pvalue)) %>%
  mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
  select(isogroup, mutated_p_updown) %>%
  na.omit()

nrow(go_input_BvsD)
head(go_input_BvsD)
nrow(go_input_BvsD)
colnames(go_input_BvsD) <- c("gene", "pval")
head(go_input_BvsD)
write.csv(go_input_BvsD, file="BvsD_GO.csv", quote=F, row.names=FALSE)
##




#################### C vs D   
colData$CvsD<-factor(colData$treat, levels=c("LFH","LFNH"))
##second term is the "control"
resCvsD <- results(dds, contrast=c("treat","LFH","LFNH"))
#how many FDR < 10%?
table(resCvsD$padj<0.1)
table(resCvsD$padj<0.05)
table(resCvsD$padj<0.01)
# 0.1= 1184
# 0.05=530
# 0.01=17
summary(resCvsD)

nrow(resCvsD[resCvsD$padj<0.05 & !is.na(resCvsD$padj),])  # Num significantly differentially expressed genes excluding the no/low count genes   #228
#530

plotMA(resCvsD, main="LFH vs LFNH")
plotMA(resCvsD, main="LFH vs LFNH + y-lim", ylim=c(-2,2))

results <- as.data.frame(resCvsD)
head(results)

nrow(resCvsD[resCvsD$padj<0.05 & resCvsD$log2FoldChange > 0 & !is.na(resCvsD$padj),])
nrow(resCvsD[resCvsD$padj<0.05 & resCvsD$log2FoldChange < 0 & !is.na(resCvsD$padj),])
#UP in C 1401
#DOWN in C 213

write.table(resCvsD, file="CvsD.txt", quote=F, sep="\t")

cd <- read.table("CvsD.txt")
#Add annotations
ID2gene<-read.csv(file = "ID2gene.tab", sep="\t", row.names=2)
ID2gene$row.names<-NULL
newcd<-transform(merge(cd,ID2gene,by=0), row.names=Row.names, Row.names=NULL)
nrow(newcd[newcd$padj<0.05 & newcd$log2FoldChange > 0 & !is.na(newcd$padj),])
nrow(newcd[newcd$padj<0.05 & newcd$log2FoldChange < 0 & !is.na(newcd$padj),])


write.table(newcd, file="CvsD.txt")




MAttPlot <- function(df) {
  df$dotcol <- ifelse(df$log2FoldChange > 0 & df$padj < 0.1, 'darkorange',
                      ifelse(df$log2FoldChange < 0 & df$padj < 0.1, 'cyan3', 'black'))
  df$baseMean <- log2(df$baseMean)
  print(head(df))
  gg <- ggplot(df, aes(baseMean, log2FoldChange)) +
    geom_point(size = .3, color = df$dotcol) +
    theme_bw() +
    theme(panel.grid = element_blank())
  print(gg)
}
MAttPlot(cd)
##########make the GO table for MWU
head(cd)
cd$isogroup=row.names(cd)
library(dplyr)
go_input_CvsD = cd %>%
  mutate(mutated_p = -log(pvalue)) %>%
  mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
  select(isogroup, mutated_p_updown) %>%
  na.omit()

nrow(go_input_CvsD)
head(go_input_CvsD)
nrow(go_input_CvsD)
colnames(go_input_CvsD) <- c("gene", "pval")
head(go_input_CvsD)
write.csv(go_input_CvsD, file="CvsD_GO.csv", quote=F, row.names=FALSE)


###############################################################################################
###############################################################################################
###############################################################################################
#--------------get pvals

#A
valAvsB=cbind(resAvsB$pvalue, resAvsB$padj)
head(valAvsB)
colnames(valAvsB)=c("HFH", "HFNH")
length(valAvsB[,1])
table(complete.cases(valAvsB))
valAvsC=cbind(resAvsC$pvalue, resAvsC$padj)
head(valAvsC)
colnames(valAvsC)=c("HFH", "LFH")
length(valAvsC[,1])
table(complete.cases(valAvsC))


#B

valBvsD=cbind(resBvsD$pvalue, resBvsD$padj)
head(valBvsD)
colnames(valBvsD)=c("HFNH", "LFNH")
length(valBvsD[,1])
table(complete.cases(valBvsD))
valBvsD
#C

valCvsD=cbind(resCvsD$pvalue, resCvsD$padj)
head(valCvsD)
colnames(valCvsD)=c("LFH", "LFNH")
length(valCvsD[,1])
table(complete.cases(valCvsD))


#################################################################################
#################################################################################
#################################################################################
####### Make rlogdata and pvals table
library(gplots)

rlog=rlogTransformation(dds, blind=TRUE) 
rld=assay(rlog)
head(rld)
colData$treatcol=paste(colData$treat, colData$colony)
colnames(rld)=paste(colData$treatcol)

head(rld)
length(rld[,1])
#HFNH vs LFNH only
rld=data.frame(rld)
colData$treatcol
HFNHvsLFNHrld= rld %>% select("HFNH.Col_a","HFNH.Col_b", "HFNH.Col_c", "HFNH.Col_d","LFNH.Col_b","LFNH.Col_a","LFNH.Col_c")

write.csv(HFNHvsLFNHrld, "HFNHvsLFNHrld.csv", quote=F)
HFNHvsLFNHrld=read.csv("HFNHvsLFNHrld.csv", row.names = 1)
head(matrix(HFNHvsLFNHrld))

#ribosome
GO=read.table(file="CC_BvsD_GO_T.csv", header=T)

ribo=filter(GO, term=="GO:0005840")
nrow(ribo[abs(ribo$value)>1,])
sigribo=ribo[abs(ribo$value)>1,]
hm=subset(HFNHvsLFNHrld, rownames(HFNHvsLFNHrld) %in% sigribo$seq)

preribo=filter(GO,term=="GO:0030684")
nrow(preribo[abs(preribo$value)>1,])
sigpreribo=preribo[abs(preribo$value)>1,]
hm=subset(HFNHvsLFNHrld, rownames(HFNHvsLFNHrld) %in% sigpreribo$seq)
#cilium and motile cilum
cilium=filter(GO,term=="GO:0005929" | term== "GO:0031514")
nrow(cilium[abs(cilium$value)>1,])
sigcilium=cilium[abs(cilium$value)>1,]
hm=subset(HFNHvsLFNHrld, rownames(HFNHvsLFNHrld) %in% sigcilium$seq)

GO=read.table(file="BP_BvsD_GO_T.csv", header=T)
#glucose/ glycolysis
#glucose metabolic process
#GO:0006006
#glycolytic process 
#GO:0006096 
glucglyc=filter(GO,term=="GO:0006006" | term== "GO:0006096")
nrow(glucglyc[abs(glucglyc$value)>1,])
sigglucglyc=glucglyc[abs(glucglyc$value)>1,]
hm=subset(HFNHvsLFNHrld, rownames(HFNHvsLFNHrld) %in% sigglucglyc$seq)

#response to oxidative stress
#GO:0006979
oxid=filter(GO,term=="GO:0006979")
nrow(oxid[abs(oxid$value)>1,])
sigoxid=oxid[abs(oxid$value)>1,]
hm=subset(HFNHvsLFNHrld, rownames(HFNHvsLFNHrld) %in% sigoxid$seq)

#carb - so many terms
#carbohydrate transport GO:0008643
#carbohydrate biosynthetic process GO:0016051
carb=filter(GO,term=="GO:0008643"| term=="GO:0016051")
nrow(carb[abs(carb$value)>1,])
sigcarb=carb[abs(carb$value)>1,]
hm=subset(HFNHvsLFNHrld, rownames(HFNHvsLFNHrld) %in% sigcarb$seq)
#phosphorylation 
#GO:0006468;GO:0016310
phos=filter(GO,term=="GO:0006468;GO:0016310")
nrow(phos[abs(phos$value)>1,])
sigphos=phos[abs(phos$value)>1,]
hm=subset(HFNHvsLFNHrld, rownames(HFNHvsLFNHrld) %in% sigphos$seq)

#mitochondrial respiratory chain complex assembly
#mitochondrial transport
#mitochondrion organization
mito=filter(GO,term=="GO:0033108;GO:0070271"| term=="GO:0006839"| term=="GO:0007005")
nrow(mito[abs(mito$value)>1,])
sigmito=mito[abs(mito$value)>1,]
hm=subset(HFNHvsLFNHrld, rownames(HFNHvsLFNHrld) %in% sigmito$seq)

#fatty acid metabolic
#GO:0006631
fatty=filter(GO,term=="GO:0006631")
nrow(fatty[abs(fatty$value)>1,])
sigfatty=fatty[abs(fatty$value)>1,]
hm=subset(HFNHvsLFNHrld, rownames(HFNHvsLFNHrld) %in% sigfatty$seq)
#tricarboxylic acid cycle
tricarb=filter(GO,term=="GO:0006099")
nrow(tricarb[abs(tricarb$value)>1,])
sigtricarb=tricarb[abs(tricarb$value)>1,]
hm=subset(HFNHvsLFNHrld, rownames(HFNHvsLFNHrld) %in% sigtricarb$seq)
#nitrogen compound transport
nitrotrans=filter(GO,term=="GO:0071705")
nrow(nitrotrans[abs(nitrotrans$value)>1,])
signitrotrans=nitrotrans[abs(nitrotrans$value)>1,]
hm=subset(HFNHvsLFNHrld, rownames(HFNHvsLFNHrld) %in% signitrotrans$seq)
#cell-matrix adhesion & extracellular matrix organization
cellmatrix=filter(GO,term=="GO:0007160;GO:0031589"| term=="GO:0030198;GO:0043062")
nrow(cellmatrix[abs(cellmatrix$value)>1,])
sigcellmatrix=cellmatrix[abs(cellmatrix$value)>1,]
hm=subset(HFNHvsLFNHrld, rownames(HFNHvsLFNHrld) %in% sigcellmatrix$seq)
#aminoacid 
#cellular amino acid biosynthetic process
aminoacid=filter(GO,term=="GO:1901607;GO:0008652")
nrow(aminoacid[abs(aminoacid$value)>1,])
sigaminoacid=aminoacid[abs(aminoacid$value)>1,]
hm=subset(HFNHvsLFNHrld, rownames(HFNHvsLFNHrld) %in% sigaminoacid$seq)
#immune response GO:0006955
immune=filter(GO,term=="GO:0006955")
nrow(immune[abs(immune$value)>1,])
sigimmune=immune[abs(immune$value)>1,]
hm=subset(HFNHvsLFNHrld, rownames(HFNHvsLFNHrld) %in% sigimmune$seq)


dev.off() 

annotations <- read.table("BvsD_T.txt")
hm1<-transform(merge(hm,annotations,by=0))
hm.z = data.matrix(hm)
hm.z = sweep(hm.z, 1L, rowMeans(hm.z), check.margin = FALSE)
hm.z.sx = apply(hm.z, 1L, sd)
hm.z = sweep(hm.z, 1L, hm.z.sx, "/", check.margin = FALSE)
hm.z = data.matrix(hm.z)
colour = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
heatmap.2(hm.z, col = colour, Rowv = TRUE, Colv = TRUE, scale = "row", 
          dendrogram = "both",
          trace = "none", labRow = hm1$gene,cexRow=.7,
          margin = c(8,38))



####
#HS comp
GO=read.table(file="BP_AvsC_GO_T.csv", header=T)
#########################################################################################
#########################################################################################
#########################################################################################
#Sample distance heatmap
library(RColorBrewer)

# Sample distance heatmap
sampleDists <- as.matrix(dist(t(rld)))
library(gplots)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          margin=c(10, 10), main="Sample Distance Matrix")


#################################################################################
#################################################################################
#################################################################################
###########################heat map of sample distances 
library(pheatmap)
library(vegan)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(dplyr)
rldpvals <- read.csv(file="Total_Host_ex_Pvals.csv", row.names=1)
head(rldpvals)
rld=rldpvals[,1:15]
head(rld)

sampleDists <- dist(t(rld))
sampleDistMatrix <- as.matrix(sampleDists)
treat=c( "HFH", "HFH", "HFH" , "HFH", "HFNH","HFNH", "HFNH", "HFNH", "LFH","LFH","LFH","LFH","LFNH","LFNH","LFNH")
colnames(sampleDistMatrix)=paste(treat)
rownames(sampleDistMatrix)=paste(treat)
#Set colors and parameters
heat.colors = colorRampPalette(rev(c("red","yellow","black")),bias=0.3)(100)
pheatmap(sampleDistMatrix,color = heat.colors,cex=0.9,border_color=NA,cluster_rows=T,cluster_cols=T)

##########################################################
##########################################################
##########################################################
#Principle Component Analysis
# rld_LF=rld[, -grep("HF", colnames(rld))]
# rld_NH=rld[, -grep("FH", colnames(rld))]
# rld_HF=rld[, -grep("LF", colnames(rld))]
# 
# rld_t=t(rld_LF)
# rld_t=t(rld_NH)
# rld_t=t(rld_HF)
rld_t=t(rld)
which(apply(rld_t, 2, var)==0) #identify which columns have zero variance, PCA will not work on these
rld_t = rld_t[ , apply(rld_t, 2, var) != 0] #removes columns with zero variance


pca <- prcomp(rld_t,center = TRUE, scale. = TRUE, na.action = na.omit) #Scale = TRUE is advisable
head(pca)
#Introduce principle components
li <- pca$sdev^2 / sum(pca$sdev^2)
pc1v <- round(li[1] * 100, 1)
pc2v <- round(li[2] * 100, 1)
pca_s <- as.data.frame(pca$x)
head(pca_s)
pca_s <- pca_s[,c(1,2)]
pca_s$Samples = row.names(pca_s)
pca_s$treat=colData$treat
#for pooled heat
pca_s$heat=colData$heat
pca_s$flow=colData$flow
pca_s$treat=colData$group

# LF_colData=subset(colData, treat=="LFH"|treat=="LFNH")
# pca_s$treat=LF_colData$treat
# NH_colData=subset(colData, treat=="HFNH"|treat=="LFNH")
# pca_s$treat=NH_colData$treat
# HF_colData=subset(colData, treat=="HFNH"|treat=="HFH")
# pca_s$treat=HF_colData$treat

head(pca_s)
library(dplyr)
library(ggplot2)
library(ggrepel)
#for pooled stuff
hull_cyl <- pca_s %>%
  group_by(heat) %>%
  dplyr::slice(chull(PC1, PC2))
colony=as.factor (c( "B", "A", "C", "D", "A","B", "C", "D", "B","A","C","D","B","A","C"))
ggplot(pca_s, aes(PC1, PC2, color = flow, pch = treat)) +
  geom_point(size=3, color="black") +
  geom_text_repel(aes(label=colony),
                  colour = "black",size=6) +
  theme_bw() +
  #stat_ellipse(geom="polygon",level=0.65, alpha = 1/3, aes(fill = heat),show.legend = NA)+
  stat_ellipse(geom="polygon",level=0.65, alpha = 1/3, aes(fill = flow),show.legend = NA)+
  xlab(paste0("PC1: ",pc1v,"% variance")) +
  ylab(paste0("PC2: ",pc2v,"% variance")) +
  #scale_fill_discrete(name = "Cond", labels = c("High Flow + Heat", "High Flow + No Heat", "Low Flow + Heat","Low Flow + No Heat"))+
  scale_colour_manual(values=cbPalette)+
  scale_fill_manual(values=cbPalette,name = "Condition", labels = c("High Flow + Heat", "High Flow + No Heat", "Low Flow + Heat","Low Flow + No Heat") )+
  scale_shape_discrete(name = "Condition", labels = c("High Flow + Heat", "High Flow + No Heat", "Low Flow + Heat","Low Flow + No Heat"))+
  ggtitle("PCoA Host",expression(paste(italic("Ex Situ"))))
#saved as us legal, landscape 
library(vegan)
adonis(pca$x ~ heat+flow+colony, data = pca_s, method='man', na.rm = TRUE)

# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
# heat       1    215403  215403 2.86079 0.17577  0.009 **
# flow       1    178248  178248 2.36734 0.14545  0.038 * 
# colony     3    154152   51384 0.68244 0.12579  0.838   
# Residuals  9    677653   75295         0.55298          
# Total     14   1225456                 1.00000          
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


#Set parameters and components
cbPalette <- c("red", "blue", "black", "purple")

hull_cyl <- pca_s %>%
  group_by(treat) %>%
  dplyr::slice(chull(PC1, PC2))
colony=as.factor (c( "B", "A", "C", "D", "A","B", "C", "D", "B","A","C","D","B","A","C"))
ggplot(pca_s, aes(PC1, PC2, color = treat)) +
  geom_point(size=3, aes(shape=colony)) +
# geom_text_repel(aes(label=colony),
 #                colour = "black",size=6) +
  theme_bw() +
  theme(text = element_text(size = 20))+
  stat_ellipse(geom="polygon",level=0.65, alpha = 1/5, aes(fill = treat),show.legend = NA)+
  xlab(paste0("PC1: ",pc1v,"% variance")) +
  ylab(paste0("PC2: ",pc2v,"% variance")) +
  #scale_fill_discrete(name = "Cond", labels = c("High Flow + Heat", "High Flow + No Heat", "Low Flow + Heat","Low Flow + No Heat"))+
  scale_colour_manual(values=cbPalette, name="Condition", labels=c("High Flow + Heat", "High Flow + No Heat", "Low Flow + Heat","Low Flow + No Heat") )+
  scale_fill_manual(values=cbPalette,name = "Condition", labels = c("High Flow + Heat", "High Flow + No Heat", "Low Flow + Heat","Low Flow + No Heat") )+
   scale_shape_manual(values=c('A' = 17, 'B' = 16, 'C'=15, 'D'=11), name = "Genotype", labels = c("A", "B", "C","D"))+
  ggtitle("PCoA Host",expression(paste(italic("Ex Situ"))))
#saved as us legal, landscape 

#Run adonis for significance
library(vegan)
adonis(pca$x ~ treat, data = pca_s, method='man', na.rm = TRUE)
library(RVAideMemoire)
?pairwise.perm.manova
#LF comp
#            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# treat      1     22554   22555  1.6791 0.25139  0.247
# Residuals  5     67164   13433         0.74861       
# Total      6     89718                 1.00000  
#NH comp
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# treat      1    131794  131794  4.3147 0.46322  0.032 *
# Residuals  5    152725   30545         0.53678         
# Total      6    284519                 1.00000         

# HF comp
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# treat      1     79416   79416  1.3784 0.18682   0.19
# Residuals  6    345688   57615         0.81318       
# Total      7    425104                 1.00000       

#all comp
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# treat      3     84119   28040  2.4631 0.40183  0.014 *
# Residuals 11    125223   11384         0.59817         
# Total     14    209342                 1.00000 

#treat sig, colony not sig


rl=rlog(dds)
e=ExpressionSet(assay(rl), AnnotatedDataFrame(as.data.frame(colData(rl))))
save(dds,colData,m,rl,file="initial.RData")
# generating normalized variance-stabilized data for PCoA and heatmaps
load("initial.RData")
vsd=assay(rl)
snames=paste(colnames(m),g,sep=".")
colnames(vsd)=snames
save(vsd,g,file="vsd.RData")
# heatmap and hierarchical clustering:
load("vsd.RData")
# similarity among samples
pheatmap(cor(vsd))

# Principle coordiante analysis
conditions=g
g=paste(treat,sep=".")
#colors and subsetters:
indcols=rep("black",ncol(vsd))
indcols[g=="HFH"]="purple"
indcols[g=="HFNH"]="red"
indcols[g=="LFH"]="green2"
indcols[g=="LFNH"]="blue"
dds.pcoa=capscale(dist(t(vsd),method="manhattan")~1)

# how many good PC's we have? 
plot(dds.pcoa$CA$eig/sum(dds.pcoa$CA$eig)) 
scores=dds.pcoa$CA$u

scores
write.csv(scores, file="ex_situ_host_scores.csv")


#####PCoA in ggplot r########
treat=c( "HFH", "HFH", "HFH" , "HFH", "HFNH","HFNH", "HFNH", "HFNH", "LFH","LFH","LFH","LFH","LFNH","LFNH","LFNH")
colony=as.factor (c( "Col_b", "Col_a", "Col_c", "Col_d", "Col_a","Col_b", "Col_c", "Col_d", "Col_b","Col_a","Col_c","Col_d","Col_b","Col_a","Col_c"))
a.con=data.frame(cbind(colony,treat))
a.con
vsd=assay(rl)
snames=paste(colnames(m),a.con[,2],sep=".")
colnames(txi.rsem[["counts"]])
colnames(vsd)=snames
save(vsd,a.con,file="vsd.RData")
load("vsd.RData")
conditions=a.con
conditions$ibyt=paste(conditions[,1],conditions[,2],sep=".")
dds.pcoa=capscale(dist(t(vsd),method="manhattan")~1)
scores=dds.pcoa$CA$u

data.scores <- as.data.frame(scores)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$sample <- rownames(data.scores)  # create a column of sample names, from the rownames of data.scores
data.scores$grp <- a.con$treat #  add group variable
data.scores$colony <- a.con$colony
library(dplyr)
hull_cyl <- data.scores %>%
  group_by(grp) %>%
  dplyr::slice(chull(MDS1, MDS2))

p <- ggplot(data.scores, aes(MDS1, MDS2)) + geom_point(shape = 21)
p + aes(fill = factor(grp)) + geom_polygon(data = hull_cyl, alpha = 0.5, color='black')+
  geom_text(data=data.scores,
            aes(MDS1, MDS2 , label = colony),
            color= "black", size=7, fontface = 2,
            hjust = 0.5, vjust = 0.5)+
  scale_fill_discrete(name = "Cond", labels = c("High Flow + Heat", "High Flow + No Heat", "Low Flow + Heat","Low Flow + No Heat"))+
  labs(title="PCoA Host")

