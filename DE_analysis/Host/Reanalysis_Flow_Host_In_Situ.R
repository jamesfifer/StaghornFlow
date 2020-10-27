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
 dir <- 'C:/Users/james/Documents/Documents/JAMES Broken Computer/Guam/Raymundo/Reanalysis/In_Situ'
 
  list.files(dir)
 #samples.txt should at minimum have the directory or file name headers (e.g. HFH_rep1)
 #ALL (use this for pcoa)
  samples <- read.table(file.path(dir, "samples.txt"), header = FALSE)
  
  
  #ALL MINUS BLEACHED
  samples <- read.table(file.path(dir, "samples_minus_bleached.txt"), header = FALSE)
 #Bleached
samples <- read.table(file.path(dir, "Bleachedsamples.txt"), header = FALSE)

  
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
 #treat=c( "BFS", "BFS", "BFS" , "BNS", "BNS","BNS", "FSO", "FSO", "FSO","FSI","FSI","FSI","NSO","NSO","NSO","NSI","NSI","NSI")
 #colony=as.factor (c( "Col_1", "Col_2", "Col_3", "Col_4", "Col_5","Col_6", "Col_1", "Col_2", "Col_3","Col_1","Col_2","Col_3","Col_4","Col_5","Col_6","Col_4","Col_5","Col_6"))
 #All (use this for pcoa)
 treat=c( "FS", "FS", "FS" , "NS", "NS","NS", "FS", "FS", "FS","FS","FS","FS","NS","NS","NS","NS","NS","NS")
  colony=as.factor (c( "Col_1", "Col_2", "Col_3", "Col_4", "Col_5","Col_6", "Col_1", "Col_2", "Col_3","Col_1","Col_2","Col_3","Col_4","Col_5","Col_6","Col_4","Col_5","Col_6"))
 O.I=c("O", "O", "O","O", "O", "O","O", "O", "O","I","I","I","O","O","O","I","I","I")
 bleach.state= c( "B", "B", "B" , "B", "B","B", "U", "U", "U","U","U","U","U","U","U","U","U","U")
 g=data.frame(treat,colony,O.I, bleach.state)
 g
 colData<- g
 
 
 # 
 #All MINUS BLEACHED #use for heatmaps
 treat=c(  "FS", "FS", "FS","FS","FS","FS","NS","NS","NS","NS","NS","NS")
 colony=as.factor (c( "Col_1", "Col_2", "Col_3","Col_1","Col_2","Col_3","Col_4","Col_5","Col_6","Col_4","Col_5","Col_6"))
 O.I=c("O", "O", "O","I","I","I","O","O","O","I","I","I")
 
#Bleached only
 treat=c("FS","FS","FS","NS","NS","NS")
 colony=as.factor (c("Col_1","Col_2","Col_3","Col_4","Col_5","Col_6"))
 
 
 #Bleached
 #treat=c( "BFS", "BFS", "BFS" , "FSO", "FSO", "FSO","FSI","FSI","FSI")
 #colony=as.factor (c( "Col_1", "Col_2", "Col_3", "Col_1", "Col_2", "Col_3","Col_1", "Col_2", "Col_3"))
 #treat=c( "BNS", "BNS", "BNS" , "NSO", "NSO", "NSO","NSI","NSI","NSI")
 #colony=as.factor (c( "Col_1", "Col_2", "Col_3", "Col_1", "Col_2", "Col_3","Col_1", "Col_2", "Col_3"))
 
  g=data.frame(treat,colony)
  g
 colData<- g
 # 
 g=data.frame(treat,colony,O.I)
  g
  colData<- g
 
 #I can't use + genotype if the genotype doesn't exist in all comparisons. 
 
 head(txi.rsem$counts)
 length(txi.rsem$counts[,1])
 #27807 host
 #22495 sym
 
 #would need to define m as an integer, not sure how kosher it is to just integer these expected counts(txi.rsem$counts)
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


#Run with bleaching/unbleached samples to look at bleach state effects
#dds<-DESeqDataSetFromTximport(txi.rsem, colData, design=~treat+bleach.state)
#Run without bleaching samples to look at FSvNS effects #I think this makes more sense to not control for OI?...
  dds<-DESeqDataSetFromTximport(txi.rsem, colData, design=~treat+ O.I)
  #tested interaction
  dds<-DESeqDataSetFromTximport(txi.rsem, colData, design=~treat+ O.I+treat:O.I)
  
  
#Run with everything for PCA
  dds<-DESeqDataSetFromTximport(txi.rsem, colData, design=~treat+ O.I+bleach.state)
  
 #also ran with the following to check for interaction but adonis was insigignificant and only 2 genes were DGE for interaction
   #dds<-DESeqDataSetFromTximport(txi.rsem, colData, design=~treat+ O.I+bleach.state+treat:bleach.state)
  #Run with bleaching samples only
  dds<-DESeqDataSetFromTximport(txi.rsem, colData, design=~treat)
#error in check full rank, so add 4th column- won't be able to control for genotype in FS vs NS comps though. 
#However this does allow us to compare O and I within and accross sites (FSvNS) while controlling for differences in individual colonies




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

min(totalCounts) #220117
max(totalCounts)  # 3707369


#run deseq
dds<-DESeq(dds)

head(dds)
res<- results(dds)


#Look at dispersions plot
plotDispEsts(dds, main="Dispersion plot ")
colData$treat
colData

#interaction 1
#run with model dds<-DESeqDataSetFromTximport(txi.rsem, colData, design=~treat+ O.I+treat:O.I)
#log2(O/I(NS)) - log2(O/I(FS))
res<-results(dds,name="treatNS.O.IO")
table(res$padj<0.1)
#FALSE  TRUE 
#21723     5
table(res$padj<0.05)
# TRUE 
#   5
table(res$padj<0.01)
#4

write.table(res, file="treatNS.O.IO.txt", quote=F, sep="\t")
cd <- read.table("treatNS.O.IO.txt")
#Add annotations
ID2gene<-read.csv(file = "ID2gene.tab", sep="\t", row.names=2)
ID2gene$row.names<-NULL
#newcd<-merge(cd, ID2gene)

newcd<-transform(merge(cd,ID2gene,by=0,all=T, row.names=Row.names, Row.names=NULL))

write.table(newcd, file="treatNS.O.IO.txt")

#interaction 2
#run with model dds<-DESeqDataSetFromTximport(txi.rsem, colData, design=~treat+ O.I+bleach.state+treat:bleach.state)
res<-results(dds,name="treatNS.bleach.stateU")
#how many FDR < 10%?
table(res$padj<0.1)
# FALSE  TRUE 
# 23190     2 
table(res$padj<0.05)
# FALSE  TRUE 
# 23190     2 
table(res$padj<0.01)
# FALSE  TRUE 
# 23191     1 





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
ID2gene<-read.csv(file = "ID2gene.tab", sep="\t", row.names=2)
ID2gene$row.names<-NULL
#newcd<-merge(cd, ID2gene)

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
#run on all minus bleached
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

nrow(resAvsC[resAvsC$padj<0.05 & !is.na(resAvsC$padj),])  # Num significantly differentially expressed genes excluding the no/low count genes   #228

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
#639
#878



write.table(resAvsC, file="FSvsNS.txt")
cd <- read.table("FSvsNS.txt")
#Add annotations
ID2gene<-read.csv(file = "ID2gene.tab", sep="\t", row.names=2)
ID2gene$row.names<-NULL
#newcd<-merge(cd, ID2gene)

newcd<-transform(merge(cd,ID2gene,by=0), row.names=Row.names, Row.names=NULL)
write.table(newcd, file="FSvsNS.txt")
nrow(newcd[newcd$padj<0.05 & newcd$log2FoldChange > 0 & !is.na(newcd$padj),])
nrow(newcd[newcd$padj<0.05 & newcd$log2FoldChange < 0 & !is.na(newcd$padj),])
#534
#308


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


plotMA(resBvsD, main="HFNH vs LFNH")
plotMA(resBvsD, main="HFNH vs LFNH + y-lim", ylim=c(-2,2))

results <- as.data.frame(resBvsD)
head(results)

nrow(resBvsD[resBvsD$padj<0.05 & resBvsD$log2FoldChange > 0 & !is.na(resBvsD$padj),])
nrow(resBvsD[resBvsD$padj<0.05 & resBvsD$log2FoldChange < 0 & !is.na(resBvsD$padj),])
#20
#42

write.table(resBvsD, file="FSvsNS_B.txt")

cd <- read.table("FSvsNS_B.txt")
#Add annotations
ID2gene<-read.csv(file = "ID2gene.tab", sep="\t", row.names=2)
ID2gene$row.names<-NULL
#newcd<-merge(cd, ID2gene)

newcd<-transform(merge(cd,ID2gene,by=0), row.names=Row.names, Row.names=NULL)
nrow(newcd[newcd$padj<0.05 & newcd$log2FoldChange > 0 & !is.na(newcd$padj),])
nrow(newcd[newcd$padj<0.05 & newcd$log2FoldChange < 0 & !is.na(newcd$padj),])
#12
#26

write.table(newcd, file="FSvsNS_B.txt")



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
####### Make rlogdata 
rlog=rlogTransformation(dds, blind=TRUE)
rld=assay(rlog)
head(rld)
colData$treatcolO.I=paste(colData$treat, colData$colony, colData$O.I)
colnames(rld)=paste(colData$treatcolO.I)

head(rld)
length(rld[,1])
rld=data.frame(rld)
#FSvsNS
FSvsNS1=rld
write.csv(FSvsNS1, "FSvsNSrld.csv", quote=F)
FSvsNS=read.csv("FSvsNSrld.csv", row.names = 1)

GO=read.table(file="CC_FSvsNS_GO.csv", header=T)
#preribosome OK
#GO:0030686 GO:0030684
preribo=filter(GO, term=="GO:0030686"| term=="GO:0030684")
nrow(preribo[abs(preribo$value)>1,])
sigpreribo=preribo[abs(preribo$value)>1,]
hm=subset(FSvsNS, rownames(FSvsNS) %in% sigpreribo$seq)
#motile cilium poor
cilium=filter(GO, term=="GO:0031514")
nrow(cilium[abs(cilium$value)>1,])
sigcilium=cilium[abs(cilium$value)>1,]
hm=subset(FSvsNS, rownames(FSvsNS) %in% sigcilium$seq)

GO=read.table(file="BP_FSvsNS_GO.csv", header=T)
#tri
#clusters well 
tricarb=filter(GO, term=="GO:0006099")
nrow(tricarb[abs(tricarb$value)>1,])
sigtricarb=tricarb[abs(tricarb$value)>1,]
hm=subset(FSvsNS, rownames(FSvsNS) %in% sigtricarb$seq)

#mitochondrial transport mitochondrial organization
#clusters well 
mito=filter(GO, term=="GO:0006839"| term=="GO:0007005")
nrow(mito[abs(mito$value)>1,])
sigmito=mito[abs(mito$value)>1,]
hm=subset(FSvsNS, rownames(FSvsNS) %in% sigmito$seq)

#glucose metabolism glycolytic process meh clustering
gluc=filter(GO, term=="GO:0006096"| term=="GO:0006006")
nrow(gluc[abs(gluc$value)>1,])
siggluc=gluc[abs(gluc$value)>1,]
hm=subset(FSvsNS, rownames(FSvsNS) %in% siggluc$seq)

#phosphorylation protein phosphorylation meh
phos=filter(GO, term=="GO:0016310"| term=="GO:0006468")
nrow(phos[abs(phos$value)>1,])
sigphos=phos[abs(phos$value)>1,]
hm=subset(FSvsNS, rownames(FSvsNS) %in% sigphos$seq)

#fatty acid metabolic meh

fatty=filter(GO, term=="GO:0006631")
nrow(fatty[abs(fatty$value)>1,])
sigfatty=fatty[abs(fatty$value)>1,]
hm=subset(FSvsNS, rownames(FSvsNS) %in% sigfatty$seq)


#nitrogen compound transport meh
nitrotrans=filter(GO, term=="GO:0071705")
nrow(nitrotrans[abs(nitrotrans$value)>1,])
signitrotrans=nitrotrans[abs(nitrotrans$value)>1,]
hm=subset(FSvsNS, rownames(FSvsNS) %in% signitrotrans$seq)

#cell-matrix adhesion extracellular matrix organization bad
cellmat=filter(GO, term=="GO:0007160;GO:0031589"| term=="GO:0030198;GO:0043062")
nrow(cellmat[abs(cellmat$value)>1,])
sigcellmat=cellmat[abs(cellmat$value)>1,]
hm=subset(FSvsNS, rownames(FSvsNS) %in% sigcellmat$seq)

#cell-cell adhesion OK
celladhesion=filter(GO, term=="GO:0007156;GO:0098609")
nrow(celladhesion[abs(celladhesion$value)>1,])
sigcelladhesion=celladhesion[abs(celladhesion$value)>1,]
hm=subset(FSvsNS, rownames(FSvsNS) %in% sigcelladhesion$seq)

#amino acid biosynthesis OK (one FS on other side)
amino=filter(GO, term=="GO:1901607;GO:0008652")
nrow(amino[abs(amino$value)>1,])
sigamino=amino[abs(amino$value)>1,]
hm=subset(FSvsNS, rownames(FSvsNS) %in% sigamino$seq)

#oxidative stress OK
oxid=filter(GO, term=="GO:0006979")
nrow(oxid[abs(oxid$value)>1,])
sigoxid=oxid[abs(oxid$value)>1,]
hm=subset(FSvsNS, rownames(FSvsNS) %in% sigoxid$seq)

#immune response good!
immune=filter(GO, term=="GO:0006955")
nrow(immune[abs(immune$value)>1,])
sigimmune=immune[abs(immune$value)>1,]
hm=subset(FSvsNS, rownames(FSvsNS) %in% sigimmune$seq)

# cilium movement and cilium or flagellum-dependent cell motility
#4FS2NS without O.I
cilium2=filter(GO, term=="GO:0003341"| term=="GO:0001539")
nrow(cilium2[abs(cilium2$value)>1,])
sigcilium2=cilium2[abs(cilium2$value)>1,]
hm=subset(FSvsNS, rownames(FSvsNS) %in% sigcilium2$seq)





annotations <- read.table("FSvsNS.txt")
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






##for pcoas#

rlog=rlogTransformation(dds, blind=TRUE) 
rld=assay(rlog)
head(rld)
colnames(rld)=paste(colData$treat)
head(rld)
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
library(dplyr)
library(ggplot2)
library(ggrepel)
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
  ggtitle("PCoA Host",expression(paste(italic("In Situ"))))
##saved as us legal, landscape 

#Run adonis for significance
adonis(pca$x ~ treat+bleach.state+O.I, data = pca_s, method='man', na.rm = TRUE)
#               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# treat         1    317437  317437 1.97889 0.11332  0.050 *
# bleach.state  1    122964  122964 0.76655 0.04390  0.662  
# O.I           1    115078  115078 0.71739 0.04108  0.705  
# Residuals    14   2245766  160412         0.80170         
# Total        17   2801244                 1.00000         
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1



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
#adonis(pca$x ~ treat, data = pca_s, method='eu', na.rm = TRUE)
# how many good PC's we have? 
plot(dds.pcoa$CA$eig/sum(dds.pcoa$CA$eig))
scores=dds.pcoa$CA$u

scores


#####PCoA in ggplot r########
treat=c( "FS", "FS", "FS" , "NS", "NS","NS", "FS", "FS", "FS","FS","FS","FS","NS","NS","NS","NS","NS","NS")
colony=as.factor (c( "Col_1", "Col_2", "Col_3", "Col_4", "Col_5","Col_6", "Col_1", "Col_2", "Col_3","Col_1","Col_2","Col_3","Col_4","Col_5","Col_6","Col_4","Col_5","Col_6"))
O.I=c("O", "O", "O","O", "O", "O","O", "O", "O","I","I","I","O","O","O","I","I","I")
bleach.state= c( "B", "B", "B" , "B", "B","B", "U", "U", "U","U","U","U","U","U","U","U","U","U")

a.con=data.frame(cbind(colony,treat, O.I, bleach.state))
a.con
vsd=assay(rl)
colnames(vsd)=snames
save(vsd,a.con,file="vsd.RData")
load("vsd.RData")
conditions=a.con
conditions$ibyt=paste(conditions[,1],conditions[,2],sep=".")
dds.pcoa=capscale(dist(t(vsd),method="manhattan")~1)
scores=dds.pcoa$CA$u
write.csv(scores, file="in_situ_host_scores.csv")

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
            aes(MDS1, MDS2 , label = bleach.state),
            color= "black", size=7, fontface = 2,
            hjust = 0.5, vjust = 0.5)+
  scale_fill_discrete(name = "Cond", labels = c("High Flow + Heat", "High Flow + No Heat", "Low Flow + Heat","Low Flow + No Heat"))+
  labs(title="PCoA Host")

