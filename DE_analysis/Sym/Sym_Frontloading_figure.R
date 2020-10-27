setwd('C:/Users/james/Documents/Documents/JAMES Broken Computer/Guam/Raymundo/Reanalysis/Sym/In_Situ/')

#NOTE ONLY FSvsNS and FSvsNS_B did I fix the read.table to read.csv problem
#BvsU_FS<-read.table('BvsU_FS.txt')
#BvsU_NS<-read.table('BvsU_NS.txt')  
#BvsU<-read.table('BvsU.txt')
FSvsNS<-read.table('FSvsNS_nobleached.txt')

FSvsNS_B<-read.table('FSvsNS_B.txt')

#FS_OvsFS_I<-read.table('FS_OvsFS_I.txt')
#NS_OvsNS_I<-read.table('NS_OvsNS_I.txt')

##HF vs LF
nrow(FSvsNS[FSvsNS$padj<0.05 & FSvsNS$log2FoldChange < 0 & !is.na(FSvsNS$padj),])
nrow(FSvsNS[FSvsNS$padj<0.05 & FSvsNS$log2FoldChange > 0 & !is.na(FSvsNS$padj),])
# 38 up in LF
# 105 up in HF

##HF vs LF bleaching
nrow(FSvsNS_B[FSvsNS_B$padj<0.05 & FSvsNS_B$log2FoldChange < 0 & !is.na(FSvsNS_B$padj),])
nrow(FSvsNS_B[FSvsNS_B$padj<0.05 & FSvsNS_B$log2FoldChange > 0 & !is.na(FSvsNS_B$padj),])
#5 up in LF
#9 up in HF

#LF bleaching vs unbleached
nrow(BvsU_NS[BvsU_NS$padj<0.05 & BvsU_NS$log2FoldChange < 0 & !is.na(BvsU_NS$padj),])
nrow(BvsU_NS[BvsU_NS$padj<0.05 & BvsU_NS$log2FoldChange > 0 & !is.na(BvsU_NS$padj),])
#2 in each 


#should do
#a) treat in situ like ex situ (heat is baseline, bleaching is heat) and look for overlap in genes
#b) should maybe also look into what was up at baseline (heat) for ex situ Fs vs NS because if
#they are stressed but not equally stressed would expect some of these related genes should be DE

#a)
in.up.HF <-(FSvsNS[FSvsNS$padj<0.05 & FSvsNS$log2FoldChange > 0 & !is.na(FSvsNS$padj),])
#

#MIGHT HAVE TO GO BACK AND CHANGE THE MERGE SO I DON'T ELIMINATE GENES WITHOUT ANNOTATIONS 
setwd('C:/Users/james/Documents/Documents/JAMES Broken Computer/Guam/Raymundo/Reanalysis/Sym/Ex_Situ')
#A: HFH B: HFNH C: LFH D:LFNH
AvsB<-read.table('AvsB.txt',  row.names = 1)

BvsD<-read.table('BvsD.txt', row.names = 1)

CvsD<-read.table('CvsD.txt', row.names = 1)

#tried .1 but there were a fuck ton of genes, so decreased to .05
##HFNH vs LFNH
nrow(BvsD[BvsD$padj<0.05 & BvsD$log2FoldChange < 0 & !is.na(BvsD$padj),])
nrow(BvsD[BvsD$padj<0.05 & BvsD$log2FoldChange > 0 & !is.na(BvsD$padj),])
#45 up in LF for .05
#553 up in HF for .05

#What genes are up at baseline
ex.up.HFNH <- (BvsD[BvsD$padj<0.05 & BvsD$log2FoldChange > 0 & !is.na(BvsD$padj),])
ex.up.LFNH <- (BvsD[BvsD$padj<0.05 & BvsD$log2FoldChange < 0 & !is.na(BvsD$padj),])

#########sidebar...are any of these up at baseline for in situ high flow conditions?
overlap.bs=merge(in.up.HF, ex.up.HFNH, by="row.names")
#yup 16 of them. 
#what about baseline for in situ low flow conditions...?
in.up.LF= FSvsNS[FSvsNS$padj<0.05 & FSvsNS$log2FoldChange < 0 & !is.na(FSvsNS$padj),]
#Pyruvate dehydrogenase most  upregualted in both (involved in glycolysis and fatty acid biosynthesis)
overlap.bs2=merge(in.up.LF, ex.up.HFNH, by="row.names")
#none of the genes upregulated in situ LF are upregulated ex situ HF

overlap3=merge(in.up.LF, ex.up.LFNH, by="row.names")
#none of the genes upregulated under insitu LF are upregulated ex situ LF
overlap4=merge(in.up.HF, ex.up.LFNH, by="row.names")
#no overlap between up in situ HF and up ex situ LF


write.csv(file="overlap.csv", overlap.bs)
# rownames(overlaps)= overlaps$Row.names
# 
# remove_these <- c(rownames(overlaps))
# template<-read.csv("AvsB_GO.csv")
# rownames(template)=template$gene
# rows_to_remove <- which(row.names(template) %in% remove_these)
# template.minus <- template[-rows_to_remove,]
# template.minus$pval <- 0
# overlaps[ ,c(2:15)] <- list(NULL)
# overlaps$pval <-1
# colnames(overlaps)[colnames(overlaps)=="Row.names"] <- "gene"
# overlaps.final=rbind(overlaps,template.minus)
# write.csv(overlaps.final, file="overlaps_GO.csv", quote=F, row.names=FALSE)

#OK back to frontloading figure....
UPNSB

#HFH vs HFNH
nrow(AvsB[AvsB$padj<0.05 & AvsB$log2FoldChange < 0 & !is.na(AvsB$padj),])
nrow(AvsB[AvsB$padj<0.05 & AvsB$log2FoldChange > 0 & !is.na(AvsB$padj),])
#up in no heat 0
#Up in Heat 1


#First take all the genes that were upregulated in the LFH vs LFNH
nrow(CvsD[CvsD$padj<0.05 & CvsD$log2FoldChange < 0 & !is.na(CvsD$padj),])
nrow(CvsD[CvsD$padj<0.05 & CvsD$log2FoldChange > 0 & !is.na(CvsD$padj),])
#Up in no heat 0
#Up in heat 52 --> we're interested in these!

UPLFH<-CvsD[CvsD$padj<0.05 & CvsD$log2FoldChange > 0 & !is.na(CvsD$padj),]

#For it to be Barshis's version of frontloading we would need to filter out the genes that
#were also upregulated in the HFH vs HFNH (although I still think these genes are interesting
#especially if they for example have 10 counts under LFNH, 50 counts under HFNH and 100 counts under Heat-
#I think this expression pattern would be 'priming' not 'frontloading', more on that later. 

UPHFH<-AvsB[AvsB$padj<0.05 & AvsB$log2FoldChange > 0 & !is.na(AvsB$padj),]

#Actually ^ I think he didn't filter these out--> WHAT WAS THIS BASED ON?? I CAN'T TELL BASED 
#ON HIS WRITING...He uses different numbers though in fig 4 than when he talks about exclusive
#upregulated in MV. Actually in venn diagram he says 136 as exclusive to MVH vs MVcontrol, and
#in the graph he says 135- I think hes referencing this 136 number though. 
#option 1
remove_these <- c(rownames(UPHFH))
remove_these
#skip, since its already not in UPLFH, itll just confuse it
#rows_to_remove <- which(row.names(UPLFH) %in% remove_these)
#UPLFH.minus.UPHFH <- UPLFH[-rows_to_remove,]
# #Y axis should be these genes & expression ratio of HFNHvsLFNH
# #Just subsetted BvsD and took the logfold change to plot on Y 
# keep_these <- c(rownames(UPLFH.minus.UPHFH))

keep_these=c(rownames(UPLFH))
row_to_keep <- which(row.names(AvsB) %in% keep_these)
HFvLFcontrol <- AvsB[row_to_keep,]
#should I have taken DGEs only? that would definitetly clean it up...Looks like Barshis does not so nvm
HFvLFcontrol$FoldChange<- 2^(HFvLFcontrol$log2FoldChange)
#^While the above is what Barshis did in his paper (according to his methods), this should remove all genes from the right side
#which is not necessarily reflective of whats actually happening. So lets look incorporate genes that were higher in HFH vs HFNH too below

#option 2 what i did for frontloading fig
Upinheat=rbind(UPHFH,UPLFH)
keep_these <- c(rownames(Upinheat))
row_to_keep <- which(row.names(AvsB) %in% keep_these)
HFvLFcontrol <- AvsB[row_to_keep,]
HFvLFcontrol$FoldChange<- 2^(HFvLFcontrol$log2FoldChange)




#X axis should be these genes and ratio of HFHvsHFNH/LFHvsLFNH
woop<-cbind(AvsB$log2FoldChange,CvsD$log2FoldChange)
rownames(woop)<-rownames(AvsB)
colnames(woop)<- c("HFHvsHFNH","LFHvsLFNH")
woop<-data.frame(woop)
woop$HFHvsHFNHfoldchange<- 2^(woop$HFHvsHFNH)
woop$LFHvsLFNHfoldchange<- 2^(woop$LFHvsLFNH)
woop$HFHvLFHfoldchangeratio<-woop$HFHvsHFNHfoldchange/woop$LFHvsLFNHfoldchange
row_to_keep <- which(row.names(woop) %in% keep_these)
HFHvHFNH.LFHvLFNHratio<- woop[row_to_keep,]
plot(HFHvHFNH.LFHvLFNHratio$HFHvLFHfoldchangeratio, HFvLFcontrol$FoldChange, xlim=range(0:2),ylim=range(0:40))+
  abline(v=1,col=3,lty=3)+ abline(h=1, col=3,lty=3)

new<-data.frame(cbind(HFHvHFNH.LFHvLFNHratio$HFHvLFHfoldchangeratio, HFvLFcontrol$FoldChange))
rownames(new)<-rownames(HFvLFcontrol)
newcd<-transform(merge(HFHvHFNH.LFHvLFNHratio,HFvLFcontrol,by=0), row.names=Row.names, Row.names=NULL)

write.csv(file="Sym_frontloading.csv", newcd)
colnames(new)<-c("HFHvLFHfoldchangeratio","HFvLFControl")
  
P1<-ggplot(newcd,aes(HFHvLFHfoldchangeratio, FoldChange)) +
  scale_y_continuous(expand = c(0,0),limits = c(0,2.5))+
  scale_x_continuous(expand=c(0,0), limits=c(0,2))+
   annotate("rect", xmin = -Inf, xmax = 1, ymin = -Inf, ymax = 1, fill= "darkgrey")  + 
  annotate("rect", xmin = -Inf, xmax = 1, ymin = 1, ymax = Inf , fill= "lightgrey") + 
  annotate("rect", xmin = 1, xmax = Inf, ymin = 1, ymax = Inf , fill= "white")+
  annotate("rect", xmin = 1, xmax = Inf, ymin = -Inf, ymax = 1 , fill= "white")+
  geom_point(aes(text=gene))+theme_minimal()+geom_hline(yintercept=1, linetype="dashed", color = "black")+
  geom_vline(xintercept=1, linetype="dashed", color = "black") +
  ylab("HF to LF control ratio")+ xlab("HF to LF foldChange ratio")+
  theme(text = element_text(size = 18))+
  theme(plot.margin = unit(c(1,1,1,1), "cm"))
P1
library(plotly)
ggplotly(P1)


##Same thing but for in situ continued 

#There's not really any differential expression gene wise for unbleaching vs bleaching 
#when looking within ns or fs
nothing=BvsU_NS[BvsU_NS$padj<0.05 & BvsU_NS$log2FoldChange < 0 & !is.na(BvsU_NS$padj),]
potentials=BvsU_NS[BvsU_NS$padj<0.05 & BvsU_NS$log2FoldChange > 0 & !is.na(BvsU_NS$padj),]

# well unreviewed Thyrotroph embryonic factor TEF AWC38_SpisGen is actually frontloaded in both 
#in situ and ex situ 
#found a lot in day/night comparisons 
#seems to be highly responsive 