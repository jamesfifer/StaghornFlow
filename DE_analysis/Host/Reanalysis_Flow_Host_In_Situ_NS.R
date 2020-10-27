 source("http://bioconductor.org/biocLite.R")
 biocLite("DESeq2")

###conduct array quality metrics to detect and remove outliers

library(DESeq) 
library(affycoretools)
library(arrayQualityMetrics)
library(genefilter)
library(Biobase)
library(pheatmap)
library(vegan)
library(rgl)
library(ape)
library(VennDiagram)

#read in raw counts from rsem and do correction by length via tximport
 
 if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
 
 BiocManager::install("tximport")
 library(tximport)
 setwd('C:/Users/james/Documents/Documents/JAMES Broken Computer/Guam/Raymundo/Reanalysis/In_Situ')
# dir <- 'C:/Users/james/Documents/BOSTON/Ecological & evolutionary genomics/Final_Project'
 dir <- 'C:/Users/james/Documents/Documents/JAMES Broken Computer/Guam/Raymundo/Reanalysis/In_Situ'
 
  list.files(dir)
 #samples.txt should at minimum have the directory or file name headers (e.g. HFH_rep1)
 #ALL 
 # samples <- read.table(file.path(dir, "samples.txt"), header = FALSE)
  
  #ALL MINUS BLEACHED
  # samples <- read.table(file.path(dir, "samples_minus_bleached.txt"), header = FALSE)

  #NS ONLY
   samples <- read.table(file.path(dir, "NSsamples.txt"), header = FALSE)
  #NS not bleached only DON;T USE
  # samples <- read.table(file.path(dir, "NSsamples2.txt"))
  
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
 
 
 #FS colonies are referred to here as Col_ 1,2,3 and NS as Col_ 4,5,6. 
 #All
#  treat=c( "FS", "FS", "FS" , "NS", "NS","NS", "FS", "FS", "FS","FS","FS","FS","NS","NS","NS","NS","NS","NS")
#  colony=as.factor (c( "Col_1", "Col_2", "Col_3", "Col_4", "Col_5","Col_6", "Col_1", "Col_2", "Col_3","Col_1","Col_2","Col_3","Col_4","Col_5","Col_6","Col_4","Col_5","Col_6"))
#  O.I=c("O", "O", "O","O", "O", "O","O", "O", "O","I","I","I","O","O","O","I","I","I")
#  bleach.state= c( "B", "B", "B" , "B", "B","B", "U", "U", "U","U","U","U","U","U","U","U","U","U")
#  
#  #All MINUS BLEACHED
#  # treat=c(  "FS", "FS", "FS","FS","FS","FS","NS","NS","NS","NS","NS","NS")
#  # colony=as.factor (c( "Col_1", "Col_2", "Col_3","Col_1","Col_2","Col_3","Col_4","Col_5","Col_6","Col_4","Col_5","Col_6"))
#  # O.I=c("O", "O", "O","I","I","I","O","O","O","I","I","I")
#  # colony.n= c( "1", "2", "3","1","2","3","1","2","3","1","2","3")
#  
#  # with #Bleached
#  treat=c( "NS", "NS", "NS" , "NS", "NS", "NS","NS","NS","NS")
#  colony=as.factor (c( "Col_1", "Col_2", "Col_3", "Col_1", "Col_2", "Col_3","Col_1", "Col_2", "Col_3"))
# O.I=c("O","O","O","O","O","O","I","I","I")
# bleach.state=c("B","B","B","U","U","U","U","U","U")

#WITH bleached per Sarah's suggestion
treat=c( "OB", "OB", "OB" , "O", "O", "O","I","I","I")
colony=as.factor (c( "Col_1", "Col_2", "Col_3", "Col_1", "Col_2", "Col_3","Col_1", "Col_2", "Col_3"))


#  #without bleached DON'T USE
# treat=c("NS", "NS", "NS","NS","NS","NS")
# colony=as.factor (c("Col_1", "Col_2", "Col_3","Col_1", "Col_2", "Col_3"))
# O.I=c("O","O","O","I","I","I")

#treat=c( "BNS", "BNS", "BNS" , "NSO", "NSO", "NSO","NSI","NSI","NSI")
 #colony=as.factor (c( "Col_1", "Col_2", "Col_3", "Col_1", "Col_2", "Col_3","Col_1", "Col_2", "Col_3"))
# g=data.frame(treat,colony,O.I)
# g
# colData<- g
#  
#  g=data.frame(treat,colony,O.I,bleach.state)
#  g
#  colData<- g
 
 #new
 g=data.frame(treat,colony)
 colData=g
 #I can't use + genotype if the genotype doesn't exist in all comparisons. 
 
 head(txi.rsem$counts)
 length(txi.rsem$counts[,1])
 #27807 host
 #22495 sym
 
 #would need to define m as an integer, not sure how kosher it is to just integer these expected counts(txi.rsem$counts)
 # m<-txi.rsem$counts
 # storage.mode(m) = "integer"
 # conditions=data.frame(treat)
 # 
 # real=newCountDataSet(m,conditions) 
 # real=estimateSizeFactors(real)
 # plot(sort(sizeFactors(real))) 
 # 
 # cds=estimateDispersions(real,method="blind")
 # vsdBlind=varianceStabilizingTransformation(cds)
 # 
 #  arrayQualityMetrics(vsdBlind,intgroup=c("treat"), force=TRUE)
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


#Run without bleaching/unbleached samples 
# dds<-DESeqDataSetFromTximport(txi.rsem, colData, design=~O.I+colony)
# dds<-DESeqDataSetFromTximport(txi.rsem, colData, design=~bleach.state+colony)
# dds<-DESeqDataSetFromTximport(txi.rsem, colData, design=~O.I+bleach.state+colony)
#new
dds<-DESeqDataSetFromTximport(txi.rsem, colData, design=~treat+colony)


#Model to compare bleaching vs unbleaching overall was: dds<-DESeqDataSetFromTximport(txi.rsem, colData, design=~treat+O.I+bleach.state)
#error in check full rank, so add 4th column- won't be able to control for genotype in FS vs NS comps though. 
#However this does allow us to compare O and I within and accross sites (FSvNS) while controlling for differences in individual colonies

#While per Sarah's suggestion if I want to do Ovs I comparison I should not have bleach stat
#as a category (instead they should be their own 3 treats), I still think it makes sense when looking
#at bleached vs non-bleached to do (~O.I+bleach.state+colony)
#model.matrix(~O.I+bleach.state+colony)


 # dds is now ready for DESeq() see DESeq2 vignette


#deseq2 was not filtering out some low guys
keep <- rowSums(counts(dds) >= 10) >= 3
dds<-dds[keep,]



totalCounts=as.vector(colSums(m))
totalCounts

dev.off()
barplot(totalCounts, col=c("coral", "coral", "coral", "red", "red", "red", "red", "blue", "blue", "blue", "blue", "green", "green" ,"green", "yellow"), ylab="raw counts")

totalCounts
#HFH  HFH  HFH  HFH  HFNH  HFNH  HFNH  HFNH  LFH  LFH  LFH  LFH  LFNH  LFNH  LFNH 
#158870  436099   66011  530727 1019095  106846  119507  110857  148101  897250   82629   60037  315173   34892  210317 

min(totalCounts) #34892
max(totalCounts)  # 1019095


#run deseq
dds<-DESeq(dds)

# Likelihood ratio test
#dds_lrt <- DESeq(dds, test="LRT", reduced = ~ 1)



head(dds)
res<- results(dds)

#Look at dispersions plot
plotDispEsts(dds, main="Dispersion plot ")
colData$treat
############################# NS #######################################################
#################### NS B vs U comparison 
#ran with dds<-DESeqDataSetFromTximport(txi.rsem, colData, design=~treat+colony)
?results
colData$AvsB<-factor(colData$treat, levels=c("OB","O"))
##second term is the "control"
resAvsB <- results(dds, contrast=c("treat","OB","O"))
#resAvsB<- results(dds, name="treatNS.bleach.stateU")
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

plotMA(resAvsB, main="B vs U")
plotMA(resAvsB, main="B vs U + y-lim", ylim=c(-2,2))



results <- as.data.frame(resAvsB)
head(results)

nrow(resAvsB[resAvsB$padj<0.05 & resAvsB$log2FoldChange > 0 & !is.na(resAvsB$padj),])
nrow(resAvsB[resAvsB$padj<0.05 & resAvsB$log2FoldChange < 0 & !is.na(resAvsB$padj),])
#UP in A = 3
#DOWN in A = 2

write.table(resAvsB, file="BvsU_NS.txt", quote=F, sep="\t")

cd <- read.table("BvsU_NS.txt")
#Add annotations
ID2gene<-read.csv(file = "ID2gene.tab", sep="\t", row.names=2)
ID2gene$row.names<-NULL
#newcd<-merge(cd, ID2gene)

newcd<-transform(merge(cd,ID2gene,by=0,all=T))
nrow(newcd[newcd$padj<0.05 & newcd$log2FoldChange > 0 & !is.na(newcd$padj),])
nrow(newcd[newcd$padj<0.05 & newcd$log2FoldChange < 0 & !is.na(newcd$padj),])
#3
#1
write.table(newcd, file="BvsU_NS.txt")
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
write.csv(go_input_AvsB, file="BvsU_NS_GO.csv", quote=F, row.names=FALSE)

#using the treat+colony, all tree are different treats
#################### Outside vs Inside 
colData$AvsC<-factor(colData$treat, levels=c("O","I"))
##second term is the "control"
resAvsC <- results(dds, contrast=c("treat","O","I"))
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
#UP in A=2
#DOWN in A=2

write.table(resAvsC, file="NS_OvsNS_I.txt", quote=F, sep="\t")

cd <- read.table("NS_OvsNS_I.txt")

#Add annotations
ID2gene<-read.csv(file = "ID2gene.tab", sep="\t", row.names=2)
ID2gene$row.names<-NULL
#newcd<-merge(cd, ID2gene)

newcd<-transform(merge(cd,ID2gene,by=0), row.names=Row.names, Row.names=NULL)
nrow(newcd[newcd$padj<0.05 & newcd$log2FoldChange > 0 & !is.na(newcd$padj),])
nrow(newcd[newcd$padj<0.05 & newcd$log2FoldChange < 0 & !is.na(newcd$padj),])

write.table(newcd, file="NS_OvsNS_I.txt")
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
write.csv(go_input_AvsC, file="NS_OvsNS_I_GO.csv", quote=F, row.names=FALSE)


