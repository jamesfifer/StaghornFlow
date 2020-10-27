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
 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
 
# BiocManager::install("tximport")
 library(tximport)
 setwd('C:/Users/james/Documents/Documents/JAMES Broken Computer/Guam/Raymundo/Reanalysis/Sym/In_Situ/')
 dir <- 'C:/Users/james/Documents/Documents/JAMES Broken Computer/Guam/Raymundo/Reanalysis/Sym/In_Situ/'
 
  list.files(dir)
 #samples.txt should at minimum have the directory or file name headers (e.g. HFH_rep1)
 #samples1 <- read.table(file.path(dir, "samples.txt"), header = FALSE)
 #all samples without outlier (cond_E_rep3)
 samples <- read.table(file.path(dir, "samples_minus_outlier.txt"), header = FALSE)
 
 samples <- read.table(file.path(dir, "samples_minus_bleached_minus_outlier.txt"), header = FALSE)
 
 #Bleached
 #samples <- read.table(file.path(dir, "Bleachedsamples.txt"), header = FALSE)
 #FS ONLY
  samples <- read.table(file.path(dir, "FSsamples.txt"), header = FALSE)
 #NS ONLY
  samples <- read.table(file.path(dir, "NSsamples_minus_outlier.txt"), header = FALSE)
 
 
 #Use this code if each file is still within its own folder
 files <- file.path(dir, samples$V2, "RSEM.genes.results")
 files
 names(files)<-samples$V2
 all(file.exists(files))
 txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
 head(txi.rsem$counts)
 txi.rsem
 
 list=read.csv(file="../Ex_Situ/all_silva.out.uniqs", header=F)
 List=list$V1
 List
 txi.rsem$counts<- txi.rsem$counts[!row.names(txi.rsem$counts) %in% List, ]
 txi.rsem$abundance<- txi.rsem$abundance[!row.names(txi.rsem$abundance) %in% List, ]
 txi.rsem$length<- txi.rsem$length[!row.names(txi.rsem$length) %in% List, ]
 #Note: there are two suggested ways of importing estimates for use with differential gene expression (DGE) methods.
 #The first method, which we show below for edgeR and for DESeq2, is to use the gene-level estimated counts from the
 #quantification tools, and additionally to use the transcript-level abundance estimates to calculate a gene-level offset
 #that corrects for changes to the average transcript length across samples. The code examples below accomplish these steps for you,
 #keeping track of appropriate matrices and calculating these offsets. For edgeR you need to assign a matrix to y$offset, but the function
 #DESeqDataSetFromTximport takes care of creation of the offset for you. Let's call this method "original counts and offset".
 
 
 
  #FS colonies are referred to here as Col_ 1,2,3 and NS as Col_ 4,5,6. 
 #treat=c( "BFS", "BFS", "BFS" , "BNS", "BNS","BNS", "FSO", "FSO", "FSO","FSI","FSI","FSI","NSO","NSO","NSO","NSI","NSI","NSI")
 #colony=as.factor (c( "Col_1", "Col_2", "Col_3", "Col_4", "Col_5","Col_6", "Col_1", "Col_2", "Col_3","Col_1","Col_2","Col_3","Col_4","Col_5","Col_6","Col_4","Col_5","Col_6"))
 #All
 # treat=c( "FS", "FS", "FS" , "NS", "NS","NS", "FS", "FS", "FS","FS","FS","FS","NS","NS","NS","NS","NS","NS")
#  colony=as.factor (c( "Col_1", "Col_2", "Col_3", "Col_4", "Col_5","Col_6", "Col_1", "Col_2", "Col_3","Col_1","Col_2","Col_3","Col_4","Col_5","Col_6","Col_4","Col_5","Col_6"))
# O.I=c("O", "O", "O","O", "O", "O","O", "O", "O","I","I","I","O","O","O","I","I","I")
# bleach.state= c( "B", "B", "B" , "B", "B","B", "U", "U", "U","U","U","U","U","U","U","U","U","U")
# g=data.frame(treat, colony, O.I, bleach.state)
# colData<- g
 #All minus outlier
 treat=c( "FS", "FS", "FS" , "NS", "NS","NS", "FS", "FS", "FS","FS","FS","FS","NS","NS","NS","NS","NS")
 colony=as.factor (c( "Col_1", "Col_2", "Col_3", "Col_4", "Col_5","Col_6", "Col_1", "Col_2", "Col_3","Col_1","Col_2","Col_3","Col_4","Col_5","Col_4","Col_5","Col_6"))
 O.I=c("O", "O", "O","O", "O", "O","O", "O", "O","I","I","I","O","O","I","I","I")
 bleach.state= c( "B", "B", "B" , "B", "B","B", "U", "U", "U","U","U","U","U","U","U","U","U")
 g=data.frame(treat, colony, O.I, bleach.state)
 colData<- g
 
 # 
 #All MINUS BLEACHED minus outlier
 treat=c(  "FS", "FS", "FS","FS","FS","FS","NS","NS","NS","NS","NS")
 colony=as.factor (c( "Col_1", "Col_2", "Col_3","Col_1","Col_2","Col_3","Col_4","Col_5","Col_4","Col_5","Col_6"))
 O.I=c("O", "O", "O","I","I","I","O","O","I","I","I")
 g=data.frame(treat,colony,O.I)
 g
 colData<- g
 
 #Bleached only
 treat=c("FS","FS","FS","NS","NS","NS")
 colony=as.factor (c("Col_1","Col_2","Col_3","Col_4","Col_5","Col_6"))
 g=data.frame(treat, colony)
 colData<- g
 
 treat=c( "OB", "OB", "OB" , "O", "O", "O","I","I","I")
 colony=as.factor (c( "Col_1", "Col_2", "Col_3", "Col_1", "Col_2", "Col_3","Col_1", "Col_2", "Col_3"))
 
 #FS Only
 treat=c( "BFS", "BFS", "BFS" , "FSO", "FSO", "FSO","FSI","FSI","FSI")
 colony=as.factor (c( "Col_1", "Col_2", "Col_3", "Col_1", "Col_2", "Col_3","Col_1", "Col_2", "Col_3"))
 
 #Bleached
 #treat=c( "BFS", "BFS", "BFS" , "FSO", "FSO", "FSO","FSI","FSI","FSI")
 #colony=as.factor (c( "Col_1", "Col_2", "Col_3", "Col_1", "Col_2", "Col_3","Col_1", "Col_2", "Col_3"))
 #treat=c( "BNS", "BNS", "BNS" , "NSO", "NSO", "NSO","NSI","NSI","NSI")
 #colony=as.factor (c( "Col_1", "Col_2", "Col_3", "Col_1", "Col_2", "Col_3","Col_1", "Col_2", "Col_3"))
 
 
 
  g=data.frame(treat, colony)
 g
 colData<- g
 
 
 
 length(txi.rsem$counts[,1])
#35597 sym
 
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
 #use this for all samples minus outlier- use for pcoa
 dds<-DESeqDataSetFromTximport(txi.rsem, colData, design=~treat+O.I+ bleach.state)
 # dds is now ready for DESeq() see DESeq2 vignette
 ##Run without bleaching samples to look at FSvNS effects #I think this makes more sense to not control for OI?...
 dds<-DESeqDataSetFromTximport(txi.rsem, colData, design=~treat+ O.I)
 #for bleaching samples
 dds<-DESeqDataSetFromTximport(txi.rsem, colData, design=~treat)
 
#deseq2 was not filtering out some low guys
keep <- rowSums(counts(dds) >= 10) >= 3
dds<-dds[keep,]



totalCounts=as.vector(colSums(m))
totalCounts

dev.off()
barplot(totalCounts, col=c("coral", "coral", "coral", "coral", "red", "red", "red", "red", "blue", "blue", "blue", "blue", "green", "green" ,"green"), ylab="raw counts")

totalCounts

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
############################# A Comparisons  #######################################################
############################# Bleached Comparison #######################################################
#################### B vs U 
colData$AvsB<-factor(colData$bleach.state, levels=c("B","U"))
##second term is the "control"
resAvsB <- results(dds, contrast=c("bleach.state","B","U"))
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

nrow(resAvsB[resAvsB$padj<0.1 & resAvsB$log2FoldChange > 0 & !is.na(resAvsB$padj),])
nrow(resAvsB[resAvsB$padj<0.1 & resAvsB$log2FoldChange < 0 & !is.na(resAvsB$padj),])
#UP in A = 0
#DOWN in A = 0

write.table(resAvsB, file="BvsU.txt", quote=F, sep="\t")
cd <- read.table("BvsU.txt")
#Add annotations
ID2gene<-read.csv(file = "Sym_ID2Description.tab", sep="\t", row.names=1)
ID2gene$row.names<-NULL

newcd<-transform(merge(cd,ID2gene,by=0,all=T, row.names=Row.names, Row.names=NULL))
write.table(newcd, file="BvsU.txt")

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
write.csv(go_input_AvsB, file="BvsU_GO.csv", quote=F, row.names=FALSE)

#################### FSvsNS###################
#run on all minus bleached minus outlier 
#design ~ treat + O.I
colData$AvsC<-factor(colData$treat, levels=c("FS","NS"))
##second term is the "control"
resAvsC <- results(dds, contrast=c("treat","FS","NS"))
#how many FDR < 10%?
table(resAvsC$padj<0.1)
table(resAvsC$padj<0.05)
table(resAvsC$padj<0.01)
# 0.1= 1
# 0.05=0
# 0.01=0
summary(resAvsC)



plotMA(resAvsC, main="HFH vs LFH")
plotMA(resAvsC, main="HFH vs LFH + y-lim", ylim=c(-2,2))

results <- as.data.frame(resAvsC)
head(results)

nrow(resAvsC[resAvsC$padj<0.05 & resAvsC$log2FoldChange > 0 & !is.na(resAvsC$padj),])
nrow(resAvsC[resAvsC$padj<0.05 & resAvsC$log2FoldChange < 0 & !is.na(resAvsC$padj),])
#with design ~treat
#445
#912
#with design ~treat + O.I. 
#used this ^
#283
#97



write.table(resAvsC, file="FSvsNS_nobleached.txt")
cd <- read.table("FSvsNS_nobleached.txt")
# #Add annotations
 ID2gene<-read.csv(file = "Sym_ID2Description.tab", sep="\t", row.names=1, header = FALSE)
 newcd<-transform(merge(cd,ID2gene,by=0),  row.names=Row.names, Row.names=NULL)
colnames(newcd)[7] <- "gene"
write.table(newcd, file="FSvsNS_nobleached.txt")
nrow(newcd[newcd$padj<0.05 & newcd$log2FoldChange > 0 & !is.na(newcd$padj),])
nrow(newcd[newcd$padj<0.05 & newcd$log2FoldChange < 0 & !is.na(newcd$padj),])
#105
#37

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
write.csv(go_input_AvsC, file="FSvsNS_GO.csv", quote=F, row.names=FALSE)


#################### FS vs NS BLEACHING ONLY
colData$BvsD<-factor(colData$treat, levels=c("FS","NS"))
##second term is the "control"
resBvsD <- results(dds, contrast=c("treat","FS","NS"))
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
#20
#9

write.table(resBvsD, file="FSvsNS_B.txt")

cd <- read.table("FSvsNS_B.txt")
#Add annotations
ID2gene<-read.csv(file = "Sym_ID2Description.tab", sep="\t", row.names=1)
#newcd<-merge(cd, ID2gene)
newcd<-transform(merge(cd,ID2gene,by=0),  row.names=Row.names, Row.names=NULL)
colnames(newcd)[7] <- "gene"
write.table(newcd, file="FSvsNS_B.txt")

nrow(newcd[newcd$padj<0.05 & newcd$log2FoldChange > 0 & !is.na(newcd$padj),])
nrow(newcd[newcd$padj<0.05 & newcd$log2FoldChange < 0 & !is.na(newcd$padj),])
#9
#5
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
write.csv(go_input_BvsD, file="FSvsNS_B_GO.csv", quote=F, row.names=FALSE)


###############################################################################################
###############################################################################################
###############################################################################################
GO=read.table(file="CC_FSvsNS.no_bleached_GO.txt", header=T)
#pre ribosome up, ribosome down 
#respiratory 



####### Make rlogdata and pvals table
####### Make rlogdata for heatmaps
library(gplots)
ID2gene<-read.csv(file = "Sym_ID2Description.tab", sep="\t", row.names=1, header = FALSE)
rlog=rlogTransformation(dds, blind=TRUE)
rld=assay(rlog)
head(rld)
colData$treatcolO.I=paste(colData$treat, colData$colony, colData$O.I)
colnames(rld)=paste(colData$treatcolO.I)

head(rld)
length(rld[,1])
rld=data.frame(rld)

FSvsNS1=rld
write.csv(FSvsNS1, "FSvsNSrldSYM.csv", quote=F)
FSvsNS=read.csv("FSvsNSrldSYM.csv", row.names = 1)
GO=read.table(file="CC_FSvsNS.no_bleached_GO.txt", header=T)
GO=read.table(file="BP_FSvsNS.no_bleached_GO.txt", header=T)
#response to oxidative stress
#GO:0006979

oxidative=filter(GO, term=="GO:0006979")
nrow(oxidative[abs(oxidative$value)>1,])
sigoxidative=oxidative[abs(oxidative$value)>1,]
hm=subset(FSvsNS, rownames(FSvsNS) %in% sigoxidative$seq)


#Ribosome biogenesis GO:0042254
ribosomebiog=filter(GO, term=="GO:0042254")
nrow(ribosomebiog[abs(ribosomebiog$value)>1,])
sigribosomebiog=ribosomebiog[abs(ribosomebiog$value)>1,]
hm=subset(FSvsNS, rownames(FSvsNS) %in% sigribosomebiog$seq)

#photorespiration GO:0009853 nothing here
ribosomebiog=filter(GO, term=="GO:0009853")
nrow(ribosomebiog[abs(ribosomebiog$value)>1,])
sigribosomebiog=ribosomebiog[abs(ribosomebiog$value)>1,]
hm=subset(FSvsNS, rownames(FSvsNS) %in% sigribosomebiog$seq)

#photosynthesis OK
#GO:0015979
photosyn=filter(GO, term=="GO:0015979")
nrow(photosyn[abs(photosyn$value)>1,])
sigphotosyn=photosyn[abs(photosyn$value)>1,]
hm=subset(FSvsNS, rownames(FSvsNS) %in% sigphotosyn$seq)

#carbohydrate transport great!
#GO:0008643
carbtrans=filter(GO, term=="GO:0008643")
nrow(carbtrans[abs(carbtrans$value)>1,])
sigcarbtrans=carbtrans[abs(carbtrans$value)>1,]
hm=subset(FSvsNS, rownames(FSvsNS) %in% sigcarbtrans$seq)

#nitrogen compound transport good
#GO:0071705
nitrogentrans=filter(GO, term=="GO:0071705")
nrow(nitrogentrans[abs(nitrogentrans$value)>1,])
signitrogentrans=nitrogentrans[abs(nitrogentrans$value)>1,]
hm=subset(FSvsNS, rownames(FSvsNS) %in% signitrogentrans$seq)

#^switch to ammonium & nitrate/ nitrite transport 
nittrans=ID2gene %>%
  rownames_to_column('gene') %>%
  filter(grepl("Ammonium transporter| Nitrate transporter| Nitrite transporter",V2))%>%
  column_to_rownames('gene')
hmGO=subset(GO, GO$seq %in% rownames(nittrans))
signittrans=(hmGO[abs(hmGO$value)>1,])
hm=subset(FSvsNS, rownames(FSvsNS) %in% signittrans$seq)

#nitreductase genes OK
nitreductase=ID2gene %>%
  rownames_to_column('gene') %>%
  filter(grepl("Nitrate reductase | Nitrite reductase",V2))%>%
  column_to_rownames('gene')
hmGO=subset(GO, GO$seq %in% rownames(nitreductase))
signitreductase=(hmGO[abs(hmGO$value)>1,])
hm=subset(FSvsNS, rownames(FSvsNS) %in% signitreductase$seq)


nrow(nitrogentrans[abs(nitrogentrans$value)>1,])
signitrogentrans=nitrogentrans[abs(nitrogentrans$value)>1,]
hm=subset(FSvsNS, rownames(FSvsNS) %in% signitrogentrans$seq)

#ABC genes OK
ABC=ID2gene %>%
  rownames_to_column('gene') %>%
  filter(grepl("ABC transporter",V2))%>%
  column_to_rownames('gene')
hmGO=subset(GO, GO$seq %in% rownames(ABC))
sigABC=(hmGO[abs(hmGO$value)>1,])
hm=subset(FSvsNS, rownames(FSvsNS) %in% sigABC$seq)

#stress genes
#heat shock 90 shows opposite pattern & 110
GO:0006457

stress=ID2gene %>%
  rownames_to_column('gene') %>%
  filter(grepl("Heat shock protein 90 | Heat shock protein 110",V2))%>%
  column_to_rownames('gene')
hmGO=subset(GO, GO$seq %in% rownames(stress))
sigstress=(hmGO[abs(hmGO$value)>1,])
hm=subset(FSvsNS, rownames(FSvsNS) %in% sigstress$seq)

#ADD A RIBOSOMAL ONE 


annotations <- read.table("FSvsNS_nobleached.txt")
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


##############
#for pca

rlog=rlogTransformation(dds, blind=TRUE) 
rld=assay(rlog)
#head(rld)
colnames(rld)=paste(colData$treat)
#head(rld)
length(rld[,1])

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
pca_s$O.I=colData$O.I
pca_s$bleach.state=colData$bleach.state
pca_s$colony=colData$colony
head(pca_s)
#Set parameters and components
cbPalette <- c("red", "blue", "black", "purple")
hull_cyl <- pca_s %>%
  group_by(treat) %>%
  dplyr::slice(chull(PC1, PC2))

ggplot(pca_s, aes(PC1, PC2, color = treat)) +
  geom_point(size=3, aes(shape=bleach.state)) + #(size=3, aes(color=treat)
  theme_bw() +
  theme(text = element_text(size = 20))+
  stat_ellipse(geom="polygon",level=0.90, alpha = 1/6, aes(fill = treat),show.legend = NA)+
  xlab(paste0("PC1: ",pc1v,"% variance")) +
  ylab(paste0("PC2: ",pc2v,"% variance")) +
  #scale_fill_discrete(name = "Condition", labels = c("High Flow + Heat", "High Flow + No Heat"))+
  scale_colour_manual(values=cbPalette, name="Condition", labels=c("High Flow Site", "Low Flow Site"))+
  scale_fill_manual(values=cbPalette,name = "Condition", labels = c("High Flow Site", "Low Flow Site") )+
  scale_shape_discrete(name = "Event", labels = c("Bleaching", "Heat Stress") )+
  ggtitle("PCoA Sym",expression(paste(italic("In Situ"))))


#Run adonis for significance
adonis(pca$x ~ treat+bleach.state+O.I, data = pca_s, method='man', na.rm = TRUE)
#only treat was significant 

#               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
# treat         1    567228  567228 3.05247 0.17447  0.002 **
# bleach.state  1    182876  182876 0.98412 0.05625  0.448   
# O.I           1     85342   85342 0.45926 0.02625  0.990   
# Residuals    13   2415741  185826         0.74303          
# Total        16   3251188                 1.00000   

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
write.csv(scores, file="sym_insitu_scores.csv")
scores


#####PCoA in ggplot r########
treat=c( "FS", "FS", "FS" , "NS", "NS","NS", "FS", "FS", "FS","FS","FS","FS","NS","NS","NS","NS","NS")
colony=as.factor (c( "Col_1", "Col_2", "Col_3", "Col_4", "Col_5","Col_6", "Col_1", "Col_2", "Col_3","Col_1","Col_2","Col_3","Col_4","Col_5","Col_4","Col_5","Col_6"))
O.I=c("O", "O", "O","O", "O", "O","O", "O", "O","I","I","I","O","O","I","I","I")
bleach.state= c( "B", "B", "B" , "B", "B","B", "U", "U", "U","U","U","U","U","U","U","U","U")
a.con=data.frame(cbind(colony,treat, O.I, bleach.state))


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
  labs(title="PCoA Sym In Situ")



