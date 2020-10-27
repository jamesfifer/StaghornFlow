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
 #Bleached
samples <- read.table(file.path(dir, "Bleachedsamples.txt"), header = FALSE)
  #FS ONLY
  samples <- read.table(file.path(dir, "FSsamples.txt"), header = FALSE)
 
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
 # treat=c( "FS", "FS", "FS" , "NS", "NS","NS", "FS", "FS", "FS","FS","FS","FS","NS","NS","NS","NS","NS","NS")
 # colony=as.factor (c( "Col_1", "Col_2", "Col_3", "Col_4", "Col_5","Col_6", "Col_1", "Col_2", "Col_3","Col_1","Col_2","Col_3","Col_4","Col_5","Col_6","Col_4","Col_5","Col_6"))
 # O.I=c("O", "O", "O","O", "O", "O","O", "O", "O","I","I","I","O","O","O","I","I","I")
 # bleach.state= c( "B", "B", "B" , "B", "B","B", "U", "U", "U","U","U","U","U","U","U","U","U","U")
 # 
 #All MINUS BLEACHED
 # treat=c(  "FS", "FS", "FS","FS","FS","FS","NS","NS","NS","NS","NS","NS")
 # colony=as.factor (c( "Col_1", "Col_2", "Col_3","Col_1","Col_2","Col_3","Col_4","Col_5","Col_6","Col_4","Col_5","Col_6"))
 # O.I=c("O", "O", "O","I","I","I","O","O","O","I","I","I")
 # colony.n= c( "1", "2", "3","1","2","3","1","2","3","1","2","3")
 
 #All FS
 #per sarah's rec I'm doing three different treatments in O.I. 
 #treat=c( "FS", "FS", "FS" , "FS", "FS", "FS","FS","FS","FS")
 colony=as.factor (c( "Col_1", "Col_2", "Col_3", "Col_1", "Col_2", "Col_3","Col_1", "Col_2", "Col_3"))
 O.I=c("OB","OB","OB","O","O","O","I","I","I")
#bleach.state=c("B","B","B","U","U","U","U","U","U")
 g=data.frame(colony,O.I)
 g
 colData<- g
 
#FS no bleach
treat=c( "FS", "FS", "FS" , "FS", "FS", "FS","FS","FS","FS")
colony=as.factor (c( "Col_1", "Col_2", "Col_3", "Col_1", "Col_2", "Col_3","Col_1", "Col_2", "Col_3"))
O.I=c("O","O","O","O","O","O","I","I","I")
bleach.state=c("B","B","B","U","U","U","U","U","U")

#treat=c( "BNS", "BNS", "BNS" , "NSO", "NSO", "NSO","NSI","NSI","NSI")
 #colony=as.factor (c( "Col_1", "Col_2", "Col_3", "Col_1", "Col_2", "Col_3","Col_1", "Col_2", "Col_3"))
 
 
 g=data.frame(colony,O.I,bleach.state)
 g
 colData<- g
 
 #I can't use + genotype if the genotype doesn't exist in all comparisons. 
 
 #head(txi.rsem$counts)
 #length(txi.rsem$counts[,1])
 #27807 host
 #22495 sym
 
 #would need to define m as an integer, not sure how kosher it is to just integer these expected counts(txi.rsem$counts)
 #m<-txi.rsem$counts
 #storage.mode(m) = "integer"
 #conditions=data.frame(treat)
 
 #real=newCountDataSet(m,conditions) 
 #real=estimateSizeFactors(real)
 #plot(sort(sizeFactors(real))) 
 
 #cds=estimateDispersions(real,method="blind")
 #vsdBlind=varianceStabilizingTransformation(cds)
 
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

#Run with all 9 (O.I is made up of OB, O or I)
 dds<-DESeqDataSetFromTximport(txi.rsem, colData, design=~colony+O.I)
 
#Run without bleaching/unbleached samples 
dds<-DESeqDataSetFromTximport(txi.rsem, colData, design=~colony+O.I)
#Model to compare bleaching vs unbleaching overall was: dds<-DESeqDataSetFromTximport(txi.rsem, colData, design=~treat+O.I+bleach.state)
#error in check full rank, so add 4th column- won't be able to control for genotype in FS vs NS comps though. 
#However this does allow us to compare O and I within and accross sites (FSvNS) while controlling for differences in individual colonies

model.matrix(~O.I+bleach.state+colony)


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
colData$O.I
############################# FS  #######################################################
#################### FS B vs U comparison 
#using OB,O, I in dispersions. Colony+O.I model
?results
colData$AvsB<-factor(colData$O.I, levels=c("OB","O"))
##second term is the "control"
resAvsB <- results(dds, contrast=c("O.I","OB","O"))
#resAvsB<- results(dds, name="treatNS.bleach.stateU")
#how many FDR < 10%?
table(resAvsB$padj<0.1)
table(resAvsB$padj<0.05)
table(resAvsB$padj<0.01)
# 0.1=29
# 0.05=21
# 0.01=15
summary(resAvsB)
colnames(resAvsB)
nrow(resAvsB[resAvsB$padj<0.05 & !is.na(resAvsB$padj),])  # Num significantly differentially expressed genes excluding the no/low count genes   #228

plotMA(resAvsB, main="B vs U")
plotMA(resAvsB, main="B vs U + y-lim", ylim=c(-2,2))



results <- as.data.frame(resAvsB)
head(results)

nrow(resAvsB[resAvsB$padj<0.1 & resAvsB$log2FoldChange > 0 & !is.na(resAvsB$padj),])
nrow(resAvsB[resAvsB$padj<0.1 & resAvsB$log2FoldChange < 0 & !is.na(resAvsB$padj),])
#UP in A = 21
#DOWN in A = 8

write.table(resAvsB, file="BvsU_FS.txt", quote=F, sep="\t")

cd <- read.table("BvsU_FS.txt")
#Add annotations
ID2gene<-read.csv(file = "ID2gene.tab", sep="\t", row.names=2)
ID2gene$row.names<-NULL


newcd<-transform(merge(cd,ID2gene,by=0,all=T))
nrow(newcd[newcd$padj<0.05 & newcd$log2FoldChange > 0 & !is.na(newcd$padj),])
nrow(newcd[newcd$padj<0.05 & newcd$log2FoldChange < 0 & !is.na(newcd$padj),])
write.table(newcd, file="BvsU_FS.txt")

#17
#4
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
write.csv(go_input_AvsB, file="BvsU_FS_GO.csv", quote=F, row.names=FALSE)

#################### Outside vs in controlling for bleached state (but still including)

colData$AvsC<-factor(colData$O.I, levels=c("O","I"))
##second term is the "control"
resAvsC <- results(dds, contrast=c("O.I","O","I"))
resultsNames(dds)



colData$AvsC

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
#UP in A=1
#DOWN in A=2

write.table(resAvsC, file="FS_OvsFS_I.txt", quote=F, sep="\t")

cd <- read.table("FS_OvsFS_I.txt")

#Add annotations
ID2gene<-read.csv(file = "ID2gene.tab", sep="\t", row.names=2)
ID2gene$row.names<-NULL
#newcd<-merge(cd, ID2gene)

newcd<-transform(merge(cd,ID2gene,by=0, row.names=Row.names, Row.names=NULL))

nrow(newcd[newcd$padj<0.05 & newcd$log2FoldChange > 0 & !is.na(newcd$padj),])
nrow(newcd[newcd$padj<0.05 & newcd$log2FoldChange < 0 & !is.na(newcd$padj),])
#1
#2

write.table(newcd, file="FS_OvsFS_I.txt")

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
write.csv(go_input_AvsC, file="FSOvsFSI_GO.csv", quote=F, row.names=FALSE)


#################### B vs D 
colData$BvsD<-factor(colData$treat, levels=c("BFS","FSI"))
##second term is the "control"
resBvsD <- results(dds, contrast=c("treat","BFS","FSI"))
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

nrow(resBvsD[resBvsD$padj<0.1 & resBvsD$log2FoldChange > 0 & !is.na(resBvsD$padj),])
nrow(resBvsD[resBvsD$padj<0.1 & resBvsD$log2FoldChange < 0 & !is.na(resBvsD$padj),])
#UP in B 2597
#DOWN in B 254

write.table(resBvsD, file="BvsD.txt", quote=F, sep="\t")

cd <- read.table("BvsD.txt")
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

rlog=rlogTransformation(dds, blind=TRUE) 
rld=assay(rlog)
head(rld)
colnames(rld)=paste(colData$treat)
head(rld)
length(rld[,1])

rldpvals=cbind(rld,valAvsB, valAvsC, valBvsD,valCvsD)
head(rldpvals)
dim(rldpvals)
# 
table(complete.cases(rldpvals))
# FALSE  TRUE 
# 14421  8074 

write.csv(rldpvals, "Total_Sym_Pvals.csv", quote=F)

head(matrix(rld))
hist(matrix(rld))
colnames(rld)=paste(colData$treat)
head(rld)


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
rldpvals <- read.csv(file="Total_Sym_Pvals.csv", row.names=1)
head(rldpvals)
rld=rldpvals[,1:15]
head(rld)

sampleDists <- dist(t(rld))
sampleDistMatrix <- as.matrix(sampleDists)
treat==c( "HFH", "HFH", "HFH" , "HFH", "HFNH","HFNH", "HFNH", "HFNH", "LFH","LFH","LFH","LFH","LFNH","LFNH","LFNH")
colnames(sampleDistMatrix)=paste(treat)
rownames(sampleDistMatrix)=paste(treat)
#Set colors and parameters
heat.colors = colorRampPalette(rev(c("red","yellow","black")),bias=0.3)(100)
pheatmap(sampleDistMatrix,color = heat.colors,cex=0.9,border_color=NA,cluster_rows=T,cluster_cols=T)

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
head(pca_s)
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
adonis(pca$x ~ treat, data = pca_s, method='eu', na.rm = TRUE)

#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
#treat      3     78993   26331  1.2327 0.25161  0.184
#Residuals 11    234957   21360         0.74839       
#Total     14    313950                 1.00000      

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


#####PCoA in ggplot r########
a.con=data.frame(cbind(colony,treat))
a.con
vsd=assay(rl)
snames=paste(colnames(countData),a.con[,2],sep=".")
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

