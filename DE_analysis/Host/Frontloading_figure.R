setwd('C:/Users/james/Documents/Documents/JAMES Broken Computer/Guam/Raymundo/Reanalysis/In_Situ')

#NOTE ONLY FSvsNS and FSvsNS_B did I fix the read.table to read.csv problem
BvsU_FS<-read.table('BvsU_FS.txt')
BvsU_NS<-read.table('BvsU_NS.txt')  
BvsU<-read.table('BvsU.txt')
FSvsNS<-read.table('FSvsNS.txt')
FSvsNS_B<-read.table('FSvsNS_B.txt')
FS_OvsFS_I<-read.table('FS_OvsFS_I.txt')
NS_OvsNS_I<-read.table('NS_OvsNS_I.txt')

##HF vs LF
nrow(FSvsNS[FSvsNS$padj<0.05 & FSvsNS$log2FoldChange < 0 & !is.na(FSvsNS$padj),])
nrow(FSvsNS[FSvsNS$padj<0.05 & FSvsNS$log2FoldChange > 0 & !is.na(FSvsNS$padj),])
# 878 up in LF
# 639 up in HF

##HF vs LF bleaching
nrow(FSvsNS_B[FSvsNS_B$padj<0.05 & FSvsNS_B$log2FoldChange < 0 & !is.na(FSvsNS_B$padj),])
nrow(FSvsNS_B[FSvsNS_B$padj<0.05 & FSvsNS_B$log2FoldChange > 0 & !is.na(FSvsNS_B$padj),])
#42 up in LF
#20 up in HF

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

#Working with the genes with annotations, could do the same thing with unannotated if you wanted to
setwd('C:/Users/james/Documents/Documents/JAMES Broken Computer/Guam/Raymundo/Reanalysis/Ex_Situ')
#A: HFH B: HFNH C: LFH D:LFNH
AvsB<-read.table('AvsB.txt')
BvsD<-read.table('BvsD.txt')
CvsD<-read.table('CvsD.txt')
#tried .1 but there were a ton of genes, so decreased to .05
##HFNH vs LFNH
nrow(BvsD[BvsD$padj<0.05 & BvsD$log2FoldChange < 0 & !is.na(BvsD$padj),])
nrow(BvsD[BvsD$padj<0.05 & BvsD$log2FoldChange > 0 & !is.na(BvsD$padj),])
#116 up in LF for .05
#3460 up in HF for .05

#What genes are up at baseline
ex.up.HFNH <- (BvsD[BvsD$padj<0.05 & BvsD$log2FoldChange > 0 & !is.na(BvsD$padj),])
ex.up.LFNH <- (BvsD[BvsD$padj<0.05 & BvsD$log2FoldChange < 0 & !is.na(BvsD$padj),])

#########non-frontloading sidebar...are any of these up at baseline for in situ high flow conditions?
overlap=merge(in.up.HF, ex.up.HFNH, by="row.names")
#yup 115 of them. 
#what about baseline for in situ low flow conditions...?
in.up.LF= FSvsNS[FSvsNS$padj<0.05 & FSvsNS$log2FoldChange < 0 & !is.na(FSvsNS$padj),]
overlap2=merge(in.up.LF, ex.up.HFNH, by="row.names")
#wow only 5 of these are upregulated under the low flow in situ!!
#GO enrichment doesn't show anything

overlaps=rbind(overlap2, overlap)
rownames(overlaps)= overlaps$Row.names

#file with annotations, pvals and GO terms
#also check up in LF in both and add to overlaps
overlap3=merge(in.up.LF, ex.up.LFNH, by="row.names")
#9
overlaps=rbind(overlaps, overlap3)
rownames(overlaps)= overlaps$Row.names
write.csv(file="overlap.csv", overlaps)

template=read.table("BP_BvsD_GO.csv", header=T)
#row.names(template)=template$seq
newcd<-transform(merge(overlaps,template,by.x=0, by.y=4))


############FRONTLOADING FIGURE###############################
#OK back to frontloading figure....
UPNSB

#HFH vs HFNH
nrow(AvsB[AvsB$padj<0.05 & AvsB$log2FoldChange < 0 & !is.na(AvsB$padj),])
nrow(AvsB[AvsB$padj<0.05 & AvsB$log2FoldChange > 0 & !is.na(AvsB$padj),])
#Up in no heat 145
#Up in Heat 76


#First take all the genes that were upregulated in the LFH vs LFNH
nrow(CvsD[CvsD$padj<0.05 & CvsD$log2FoldChange < 0 & !is.na(CvsD$padj),])
nrow(CvsD[CvsD$padj<0.05 & CvsD$log2FoldChange > 0 & !is.na(CvsD$padj),])
#Up in no heat 115
#Up in heat 1100 --> we're interested in these!

UPLFH<-CvsD[CvsD$padj<0.05 & CvsD$log2FoldChange > 0 & !is.na(CvsD$padj),]

#For it to be Barshis's version of frontloading we would need to filter out the genes that
#were also upregulated in the HFH vs HFNH (although I still think these genes are interesting
#especially if they for example have 10 counts under LFNH, 50 counts under HFNH and 100 counts under Heat-
#I think this expression pattern would be 'priming' not 'frontloading', more on that later. 

UPHFH<-AvsB[AvsB$padj<0.05 & AvsB$log2FoldChange > 0 & !is.na(AvsB$padj),]
nrow(UPHFH)
#Actually ^ I think he didn't filter these out--> WHAT WAS THIS BASED ON?? I CAN'T TELL BASED 
#ON HIS WRITING...He uses different numbers though in fig 4 than when he talks about exclusive
#upregulated in MV. Actually in venn diagram he says 136 as exclusive to MVH vs MVcontrol, and
#in the graph he says 135- I think hes referencing this 136 number though. 
#Conclusion: remove those that were significantly upregulated between HFH and HFNH
#option 1 DO EITHER THIS OR OPTION 2 DEPENDING ON WHAT YOU ARE INTERESTED IN 

remove_these <- c(rownames(UPHFH))
rows_to_remove <- which(row.names(UPLFH) %in% remove_these)
UPLFH.minus.UPHFH <- UPLFH[-rows_to_remove,]
# #Y axis should be these genes & expression ratio of HFNHvsLFNH
# #Just subsetted BvsD and took the logfold change to plot on Y 
 keep_these <- c(rownames(UPLFH.minus.UPHFH))
row_to_keep <- which(row.names(AvsB) %in% keep_these)
HFvLFcontrol <- AvsB[row_to_keep,]
#should I have taken DGEs only? that would definitetly clean it up...Looks like Barshis does not so nvm
HFvLFcontrol$FoldChange<- 2^(HFvLFcontrol$log2FoldChange)
#^While the above is what Barshis did in his paper (according to his methods), this should remove all genes from the right side
#which is not necessarily reflective of whats actually happening. So lets look incorporate genes that were higher in HFH vs HFNH too below

#option 2 this option was in our paper
Upinheat=rbind(UPHFH,UPLFH)
keep_these <- c(rownames(Upinheat))
row_to_keep <- which(row.names(AvsB) %in% keep_these)
HFvLFcontrol <- AvsB[row_to_keep,]
HFvLFcontrol$FoldChange<- 2^(HFvLFcontrol$log2FoldChange)



#X axis should be these genes and ratio of HFHvsHFNH/LFHvsLFNH
library(ggplot2)
poop<-cbind(AvsB$log2FoldChange,CvsD$log2FoldChange)
rownames(poop)<-rownames(AvsB)
colnames(poop)<- c("HFHvsHFNH","LFHvsLFNH")
poop<-data.frame(poop)
poop$HFHvsHFNHfoldchange<- 2^(poop$HFHvsHFNH)
poop$LFHvsLFNHfoldchange<- 2^(poop$LFHvsLFNH)
poop$HFHvLFHfoldchangeratio<-poop$HFHvsHFNHfoldchange/poop$LFHvsLFNHfoldchange
row_to_keep <- which(row.names(poop) %in% keep_these)
HFHvHFNH.LFHvLFNHratio<- poop[row_to_keep,]
plot(HFHvHFNH.LFHvLFNHratio$HFHvLFHfoldchangeratio, HFvLFcontrol$FoldChange, xlim=range(0:2),ylim=range(0:40))+
  abline(v=1,col=3,lty=3)+ abline(h=1, col=3,lty=3)

new<-data.frame(cbind(HFHvHFNH.LFHvLFNHratio$HFHvLFHfoldchangeratio, HFvLFcontrol$FoldChange))
rownames(new)<-rownames(HFvLFcontrol)
newcd<-transform(merge(HFHvHFNH.LFHvLFNHratio,HFvLFcontrol,by=0), row.names=Row.names, Row.names=NULL)
write.csv(file="Host_frontloading.csv", newcd)

colnames(new)<-c("HFHvLFHfoldchangeratio","HFvLFControl")
  
P1<-ggplot(newcd,aes(HFHvLFHfoldchangeratio, FoldChange)) +
  scale_y_continuous(expand = c(0,0),limits = c(0,40))+
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

####looking at the pooled stuff####
setwd('C:/Users/james/Documents/Documents/JAMES Broken Computer/Guam/Raymundo/Reanalysis/Ex_Situ')
interaction=read.table('heatHeat.flowLF.txt')
heat_vs_noheat=read.table('HvsNH.txt')
