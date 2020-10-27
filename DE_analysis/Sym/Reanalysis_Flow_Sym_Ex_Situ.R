if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq")

###conduct array quality metrics to detect and remove outliers

library(DESeq) 
#BiocManager::install("affycoretools")
library(affycoretools)
#BiocManager::install("arrayQualityMetrics")

library(arrayQualityMetrics)
library(genefilter)
library(Biobase)
#install.packages("pheatmap")
library(pheatmap)
library(vegan)
#install.packages("rgl")
library(rgl)
library(ape)
library(VennDiagram)

#read in raw counts from rsem and do correction by length via tximport
 
 if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
 
 BiocManager::install("tximport")
 library(tximport)
 setwd('C:/Users/james/Documents/Documents/JAMES Broken Computer/Guam/Raymundo/Reanalysis/Sym/Ex_Situ/')
 dir <- 'C:/Users/james/Documents/Documents/JAMES Broken Computer/Guam/Raymundo/Reanalysis/Sym/Ex_Situ/'
 
  list.files(dir)
 #samples.txt should at minimum have the directory or file name headers (e.g. HFH_rep1)
 samples <- read.table(file.path(dir, "samples.txt"), header = FALSE)
 
 #Use this code if each file is still within its own folder
 files <- file.path(dir, samples$V2, "RSEM.genes.results")
 files
 names(files)<-samples$V2
 all(file.exists(files))
 txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
 head(txi.rsem$counts)
 #remove silva artifacts
 list=read.csv(file="all_silva.out.uniqs", header=F)
 List=list$V1
 List
txi.rsem$counts<- txi.rsem$counts[!row.names(txi.rsem$counts) %in% List, ]
txi.rsem$abundance<- txi.rsem$abundance[!row.names(txi.rsem$abundance) %in% List, ]
txi.rsem$length<- txi.rsem$length[!row.names(txi.rsem$length) %in% List, ]
ncol(txi.rsem$counts)

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
 
 #for heat comp and adonis
 treat=c( "Heat", "Heat", "Heat" , "Heat", "Ambient","Ambient", "Ambient", "Ambient", "Heat","Heat","Heat","Heat","Ambient","Ambient","Ambient")
 flow=c( "HF", "HF", "HF" , "HF", "HF","HF", "HF", "HF", "LF","LF","LF","LF","LF","LF","LF")
 colony=as.factor (c( "Col_b", "Col_a", "Col_c", "Col_d", "Col_a","Col_b", "Col_c", "Col_d", "Col_b","Col_a","Col_c","Col_d","Col_b","Col_a","Col_c"))
 g=data.frame(heat,flow, colony)
 g
 colData<- g
 
 
 
 
 length(txi.rsem$counts[,1])
#35549 sym WITHOUT rrna
 
 #need to define m as an integer for array quality metrics (its not changing the numbers)
 m<-txi.rsem$counts
 storage.mode(m) = "integer"
 conditions=data.frame(treat)
 
 real=newCountDataSet(m,conditions) 
 real=estimateSizeFactors(real)
 plot(sort(sizeFactors(real))) 
 
 cds=estimateDispersions(real,method="blind")
 vsdBlind=varianceStabilizingTransformation(cds)
 
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


#dds<-DESeqDataSetFromTximport(txi.rsem, colData, design=m2)
 #treat+genotype is what Michael Love uses for treat+batch to "control for batch effects"
 #
 dds<-DESeqDataSetFromTximport(txi.rsem, colData, design=~treat+colony)
 # dds is now ready for DESeq() see DESeq2 vignette
 
 #This model was for looking for signature of heat stress and for the adonis to
 #look for clustering according to flow/heat
 
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
#Host?
#HFH  HFH  HFH  HFH  HFNH  HFNH  HFNH  HFNH  LFH  LFH  LFH  LFH  LFNH  LFNH  LFNH 
#158870  436099   66011  530727 1019095  106846  119507  110857  148101  897250   82629   60037  315173   34892  210317 
#Sym 
#166368  463573   68788  575760 1142919  110592  124799  117270  155072  903005   84027   57791  319484   37397  208461
min(totalCounts) #34892       #37397
max(totalCounts)  # 1019095   #1142919


#run deseq
dds<-DESeq(dds)


# Likelihood ratio test
#dds_lrt <- DESeq(dds, test="LRT", reduced = ~ 1)



head(dds)
res<- results(dds)

#Look at dispersions plot
plotDispEsts(dds, main="Dispersion plot ")
colData$treat

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
#0
#0

write.table(resHvsNH, file="HvsNH.txt", quote=F, sep="\t")

cd <- read.table("HvsNH.txt")

#Add annotations
ID2gene<-read.csv(file = "Sym_ID2Description.tab", sep="\t", header=F, row.names = 1)
newcd<-transform(merge(cd,ID2gene,by=0))

#annotated
nrow(newcd[newcd$padj<0.05 & newcd$log2FoldChange > 0 & !is.na(newcd$padj),])
nrow(newcd[newcd$padj<0.05 & newcd$log2FoldChange < 0 & !is.na(newcd$padj),])
#0
#0


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
#448
#31

write.table(resHFvsLF, file="HFvsLF.txt", quote=F, sep="\t")

cd <- read.table("HFvsLF.txt")

#Add annotations
ID2gene<-read.csv(file = "Sym_ID2Description.tab", sep="\t", header=F, row.names = 1)

newcd<-transform(merge(cd,ID2gene,by=0))

#annotated
nrow(newcd[newcd$padj<0.05 & newcd$log2FoldChange > 0 & !is.na(newcd$padj),])
nrow(newcd[newcd$padj<0.05 & newcd$log2FoldChange < 0 & !is.na(newcd$padj),])
#237
#15

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

############################# A Comparisons  #######################################################
#################### A vs B 
colData$AvsB<-factor(colData$treat, levels=c("HFH","HFNH"))
##second term is the "control"
resAvsB <- results(dds, contrast=c("treat","HFH","HFNH"))
#how many FDR < 10%?

table(resAvsB$padj<0.1)
table(resAvsB$padj<0.05)
table(resAvsB$padj<0.01)
# 0.1=1
# 0.05=1
# 0.01=1
summary(resAvsB)
colnames(resAvsB)
nrow(resAvsB[resAvsB$padj<0.05 & !is.na(resAvsB$padj),])  # Num significantly differentially expressed genes excluding the no/low count genes   #228

plotMA(resAvsB, main="HFH vs HFNH")
plotMA(resAvsB, main="HFH vs HFNH + y-lim", ylim=c(-2,2))



results <- as.data.frame(resAvsB)
head(results)

nrow(resAvsB[resAvsB$padj<0.05 & resAvsB$log2FoldChange > 0 & !is.na(resAvsB$padj),])
nrow(resAvsB[resAvsB$padj<0.05 & resAvsB$log2FoldChange < 0 & !is.na(resAvsB$padj),])
#UP in A = 1
#DOWN in A = 0



write.table(resAvsB, file="AvsB.txt", quote=F, sep="\t")

cd <- read.table("AvsB.txt")

#Add annotations
ID2gene<-read.csv(file = "Sym_ID2Description.tab", sep="\t", header=F, row.names = 1)

newcd<-transform(merge(cd,ID2gene,by=0),  row.names=Row.names, Row.names=NULL)
colnames(newcd)[7] <- "gene"
nrow(newcd[newcd$padj<0.05 & newcd$log2FoldChange > 0 & !is.na(newcd$padj),])
nrow(newcd[newcd$padj<0.05 & newcd$log2FoldChange < 0 & !is.na(newcd$padj),])
#1
#0

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
#DOWN in A=0 

write.table(resAvsC, file="AvsC.txt", quote=F, sep="\t")
#Add annotations
cd <- read.table("AvsC.txt")
ID2gene<-read.csv(file = "Sym_ID2Description.tab", sep="\t", header=F, row.names = 1)
newcd<-transform(merge(cd,ID2gene,by=0))
nrow(newcd[newcd$padj<0.05 & newcd$log2FoldChange > 0 & !is.na(newcd$padj),])
nrow(newcd[newcd$padj<0.05 & newcd$log2FoldChange < 0 & !is.na(newcd$padj),])
#0
#0


write.table(newcd, file="AvsC.txt")

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


#################### B vs D #The better flow comparison (No heated samples)--> Used this
colData$BvsD<-factor(colData$treat, levels=c("HFNH","LFNH"))
##second term is the "control"
resBvsD <- results(dds, contrast=c("treat","HFNH","LFNH"))
#how many FDR < 10%?
table(resBvsD$padj<0.1)
table(resBvsD$padj<0.05)
table(resBvsD$padj<0.01)
# 0.1= 2851
# 0.05=1513
# 0.01=294
summary(resBvsD)

nrow(resBvsD[resBvsD$padj<0.05 & !is.na(resBvsD$padj),])  # Num significantly differentially expressed genes excluding the no/low count genes   
#1513

plotMA(resBvsD, main="HFNH vs LFNH")
plotMA(resBvsD, main="HFNH vs LFNH + y-lim", ylim=c(-2,2))

results <- as.data.frame(resBvsD)
head(results)

nrow(resBvsD[resBvsD$padj<0.05 & resBvsD$log2FoldChange > 0 & !is.na(resBvsD$padj),])
nrow(resBvsD[resBvsD$padj<0.05 & resBvsD$log2FoldChange < 0 & !is.na(resBvsD$padj),])
#UP in B 1178
#DOWN in B 109

write.table(resBvsD, file="BvsD.txt", quote=F, sep="\t")

cd <- read.table("BvsD.txt")
#Add annotations
ID2gene<-read.csv(file = "Sym_ID2Description.tab", sep="\t", header=F, row.names = 1)

newcd<-transform(merge(cd,ID2gene,by=0),  row.names=Row.names, Row.names=NULL)
colnames(newcd)[7] <- "gene"
nrow(newcd[newcd$padj<0.05 & newcd$log2FoldChange > 0 & !is.na(newcd$padj),])
nrow(newcd[newcd$padj<0.05 & newcd$log2FoldChange < 0 & !is.na(newcd$padj),])
#553
#45
write.table(newcd, file="BvsD.txt")

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



plotMA(resCvsD, main="LFH vs LFNH")
plotMA(resCvsD, main="LFH vs LFNH + y-lim", ylim=c(-2,2))

results <- as.data.frame(resCvsD)
head(results)

nrow(resCvsD[resCvsD$padj<0.05 & resCvsD$log2FoldChange > 0 & !is.na(resCvsD$padj),])
nrow(resCvsD[resCvsD$padj<0.05 & resCvsD$log2FoldChange < 0 & !is.na(resCvsD$padj),])
#UP in C 121
#DOWN in C 1

write.table(resCvsD, file="CvsD.txt", quote=F, sep="\t")

cd <- read.table("CvsD.txt")
#Add annotations
ID2gene<-read.csv(file = "Sym_ID2Description.tab", sep="\t", header=F, row.names = 1)
newcd<-transform(merge(cd,ID2gene,by=0),  row.names=Row.names, Row.names=NULL)
colnames(newcd)[7] <- "gene"
nrow(newcd[newcd$padj<0.05 & newcd$log2FoldChange > 0 & !is.na(newcd$padj),])
nrow(newcd[newcd$padj<0.05 & newcd$log2FoldChange < 0 & !is.na(newcd$padj),])
#52
#0
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

write.csv(HFNHvsLFNHrld, "SYM_HFNHvsLFNHrld.csv", quote=F)
HFNHvsLFNHrld=read.csv("SYM_HFNHvsLFNHrld.csv", row.names = 1)
head(matrix(HFNHvsLFNHrld))



GO=read.table(file="CC_BvsD_GO.csv", header=T)
#preribosome is mixed, ribosome is down, photosystem is up
#photosystem GO:0009521;GO:0009523 #didn't save, photosynthesis (BP) terms tell redundant story. 
photosystem=filter(GO, term=="GO:0009521;GO:0009523")
nrow(photosystem[abs(photosystem$value)>1,])
sigphotosystem=photosystem[abs(photosystem$value)>1,]
hm=subset(HFNHvsLFNHrld, rownames(HFNHvsLFNHrld) %in% sigphotosystem$seq)
#preribosome GO:0030684;GO:0032040
preribosome=filter(GO, term=="GO:0030684;GO:0032040")
nrow(preribosome[abs(preribosome$value)>1,])
sigpreribosome=preribosome[abs(preribosome$value)>1,]
hm=subset(HFNHvsLFNHrld, rownames(HFNHvsLFNHrld) %in% sigpreribosome$seq)

GO=read.table(file="BP_BvsD_GO.csv", header=T)
#photosynthesis overall is up
#carbon utilization
#lipid translocation up
#cellular response to nitrogen compound up
#nitrogen compound transport up (lots)
#oxidative stress up
#reactive oxygen species metabolism up
#host terms up
#nutrient




ID2gene<-read.csv(file = "Sym_ID2Description.tab", sep="\t", header=F, row.names = 1)
#photorespiration GO:0009853
ribosomebiog=filter(GO, term=="GO:0009853")
nrow(ribosomebiog[abs(ribosomebiog$value)>1,])
sigribosomebiog=ribosomebiog[abs(ribosomebiog$value)>1,]
hm=subset(FSvsNS, rownames(FSvsNS) %in% sigribosomebiog$seq)

stress=ID2gene %>%
  rownames_to_column('gene') %>%
  filter(grepl("Heat shock protein",V2))%>%
  column_to_rownames('gene')
hmGO=subset(GO, GO$seq %in% rownames(stress))
sigstress=(hmGO[abs(hmGO$value)>1,])
hm=subset(HFNHvsLFNHrld, rownames(HFNHvsLFNHrld) %in% sigstress$seq)

#response to oxidative stress
#GO:0006979

oxidative=filter(GO, term=="GO:0006979")
nrow(oxidative[abs(oxidative$value)>1,])
sigoxidative=oxidative[abs(oxidative$value)>1,]
hm=subset(HFNHvsLFNHrld, rownames(HFNHvsLFNHrld) %in% sigoxidative$seq)


#Ribosome biogenesis GO:0042254
ribosomebiog=filter(GO, term=="GO:0042254")
nrow(ribosomebiog[abs(ribosomebiog$value)>1,])
sigribosomebiog=ribosomebiog[abs(ribosomebiog$value)>1,]
hm=subset(HFNHvsLFNHrld, rownames(HFNHvsLFNHrld) %in% sigribosomebiog$seq)


#photosynthesis, photosytem is downregulated in both in situ and ex situ HF- increased dependence on heterotrophy?
#GO:0015979
photosyn=filter(GO, term=="GO:0015979")
nrow(photosyn[abs(photosyn$value)>1,])
sigphotosyn=photosyn[abs(photosyn$value)>1,]
hm=subset(HFNHvsLFNHrld, rownames(HFNHvsLFNHrld) %in% sigphotosyn$seq)

#carbohydrate transport great!
#GO:0008643
carbtrans=filter(GO, term=="GO:0008643")
nrow(carbtrans[abs(carbtrans$value)>1,])
sigcarbtrans=carbtrans[abs(carbtrans$value)>1,]
hm=subset(HFNHvsLFNHrld, rownames(HFNHvsLFNHrld) %in% sigcarbtrans$seq)



#ammonium & nitrate/ nitrite transport 
nittrans=ID2gene %>%
  rownames_to_column('gene') %>%
  filter(grepl("Ammonium transporter| Nitrate transporter| Nitrite transporter",V2))%>%
  column_to_rownames('gene')
hmGO=subset(GO, GO$seq %in% rownames(nittrans))
signittrans=(hmGO[abs(hmGO$value)>1,])
hm=subset(HFNHvsLFNHrld, rownames(HFNHvsLFNHrld) %in% signittrans$seq)

#nitreductase genes OK
nitreductase=ID2gene %>%
  rownames_to_column('gene') %>%
  filter(grepl("Nitrate reductase | Nitrite reductase",V2))%>%
  column_to_rownames('gene')
hmGO=subset(GO, GO$seq %in% rownames(nitreductase))
signitreductase=(hmGO[abs(hmGO$value)>1,])
hm=subset(HFNHvsLFNHrld, rownames(HFNHvsLFNHrld) %in% signitreductase$seq)




#ABC genes OK
ABC=ID2gene %>%
  rownames_to_column('gene') %>%
  filter(grepl("ABC transporter",V2))%>%
  column_to_rownames('gene')
hmGO=subset(GO, GO$seq %in% rownames(ABC))
sigABC=(hmGO[abs(hmGO$value)>1,])
hm=subset(HFNHvsLFNHrld, rownames(HFNHvsLFNHrld) %in% sigABC$seq)

annotations <- read.table("BvsD.txt")

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
          margin = c(8,37))



### for pcoa

rlog=rlogTransformation(dds, blind=TRUE) 
rld=assay(rlog)
#head(rld)
colnames(rld)=paste(colData$treat)
#head(rld)
length(rld[,1])


head(matrix(rld))
hist(matrix(rld))
colnames(rld)=paste(colData$treat)
head(rld)

##########################################################
##########################################################
##########################################################
#Principle Component Analysis

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
#
hull_cyl <- pca_s %>%
  group_by(treat) %>%
  dplyr::slice(chull(PC1, PC2))
colony=as.factor (c( "B", "A", "C", "D", "A","B", "C", "D", "B","A","C","D","B","A","C"))
ggplot(pca_s, aes(PC1, PC2, color = treat)) +
  geom_point(size=3, aes(shape=colony)) +
   #geom_text_repel(aes(label=Samples),
  #               colour = "black",size=6) +
  theme_bw() +
  theme(text = element_text(size = 20))+
  stat_ellipse(geom="polygon",level=0.65, alpha = 1/5, aes(fill = treat),show.legend = NA)+
  xlab(paste0("PC1: ",pc1v,"% variance")) +
  ylab(paste0("PC2: ",pc2v,"% variance")) +
  #scale_fill_discrete(name = "Cond", labels = c("High Flow + Heat", "High Flow + No Heat", "Low Flow + Heat","Low Flow + No Heat"))+
  scale_colour_manual(values=cbPalette, name="Condition", labels=c("High Flow + Heat", "High Flow + No Heat", "Low Flow + Heat","Low Flow + No Heat") )+
  scale_fill_manual(values=cbPalette,name = "Condition", labels = c("High Flow + Heat", "High Flow + No Heat", "Low Flow + Heat","Low Flow + No Heat") )+
  scale_shape_manual(values=c('A' = 17, 'B' = 16, 'C'=15, 'D'=11), name = "Genotype", labels = c("A", "B", "C","D"))+
  ggtitle("PCoA Sym",expression(paste(italic("Ex Situ"))))

#flow
hull_cyl <- pca_s %>%
  group_by(treat) %>%
  dplyr::slice(chull(PC1, PC2))
colony=as.factor (c( "B", "A", "C", "D", "A","B", "C", "D", "B","A","C","D","B","A","C"))
ggplot(pca_s, aes(PC1, PC2, color = flow)) +
  geom_point(size=3, aes(shape=colony)) +
  #geom_text_repel(aes(label=Samples),
  #               colour = "black",size=6) +
  theme_bw() +
  theme(text = element_text(size = 20))+
  stat_ellipse(geom="polygon",level=0.65, alpha = 1/5, aes(fill = flow),show.legend = NA)+
  xlab(paste0("PC1: ",pc1v,"% variance")) +
  ylab(paste0("PC2: ",pc2v,"% variance")) +
  #scale_fill_discrete(name = "Cond", labels = c("High Flow + Heat", "High Flow + No Heat", "Low Flow + Heat","Low Flow + No Heat"))+
  scale_colour_manual(values=cbPalette, name="Condition", labels=c("High Flow", "Low Flow") )+
  scale_fill_manual(values=cbPalette,name = "Condition", labels = c("High Flow", "Low Flow") )+
  scale_shape_manual(values=c('A' = 17, 'B' = 16, 'C'=15, 'D'=11), name = "Genotype", labels = c("A", "B", "C","D"))+
  ggtitle("PCoA Sym",expression(paste(italic("Ex Situ"))))



#Set parameters and components
cbPalette <- c("red", "blue", "black", "purple")
ggplot(pca_s, aes(PC1, PC2, color = treat, pch = treat)) +
  geom_point(size=3) +
  geom_text_repel(aes(label=Samples)) +
  scale_colour_manual(values=cbPalette)+
  theme_bw() +
  # geom_density2d(alpha=.5)+
  geom_polygon(alpha=.2)+
  xlab(paste0("PC1: ",pc1v,"% variance")) +
  ylab(paste0("PC2: ",pc2v,"% variance")) 

#Run adonis for significance
adonis(pca$x ~ heat+flow+colony, data = pca_s, method='eu', na.rm = TRUE)

#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# heat       1      5546  5546.1 0.64634 0.04299  0.786  
# flow       1     23427 23427.4 2.73021 0.18159  0.013 *
# colony     3     22809  7603.0 0.88605 0.17680  0.608  
# Residuals  9     77227  8580.8         0.59862         
# Total     14    129010                 1.00000

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
write.csv(scores, file="ex_situ_sym_scores.csv")


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
  labs(title="PCoA Sym Ex Situ")



