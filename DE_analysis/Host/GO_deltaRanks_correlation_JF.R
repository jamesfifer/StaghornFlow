
####cleaned up function (dirty is at bottom)
#
#GO.csv has to be created first using GOMWU script 
BP_HFNHvsLFNH=read.table("~/Documents/JAMES Broken Computer/Guam/Raymundo/Reanalysis/Ex_Situ/MWU_BP_BvsD_GO.csv",header=T) #ambient temp comp
MF_HFNHvsLFNH=read.table("~/Documents/JAMES Broken Computer/Guam/Raymundo/Reanalysis/Ex_Situ/MWU_MF_BvsD_GO.csv",header=T) #ambient temp comp
CC_HFNHvsLFNH=read.table("~/Documents/JAMES Broken Computer/Guam/Raymundo/Reanalysis/Ex_Situ/MWU_CC_BvsD_GO.csv",header=T) #ambient temp comp
HFNHvsLFNH=rbind(BP_HFNHvsLFNH,MF_HFNHvsLFNH,CC_HFNHvsLFNH)

BP_HFHvsLFH=read.table("~/Documents/JAMES Broken Computer/Guam/Raymundo/Reanalysis/Ex_Situ/MWU_BP_AvsC_GO.csv",header=T)
MF_HFHvsLFH=read.table("~/Documents/JAMES Broken Computer/Guam/Raymundo/Reanalysis/Ex_Situ/MWU_MF_AvsC_GO.csv",header=T)
CC_HFHvsLFH=read.table("~/Documents/JAMES Broken Computer/Guam/Raymundo/Reanalysis/Ex_Situ/MWU_CC_AvsC_GO.csv",header=T)
HFHvsLFH=rbind(BP_HFHvsLFH,MF_HFHvsLFH,CC_HFHvsLFH)

BP_NSOvsNSI=read.table("~/Documents/JAMES Broken Computer/Guam/Raymundo/Reanalysis/In_Situ/MWU_BP_NS_OvsNS_I_GO.csv", header=T)
MF_NSOvsNSI=read.table("~/Documents/JAMES Broken Computer/Guam/Raymundo/Reanalysis/In_Situ/MWU_MF_NS_OvsNS_I_GO.csv", header=T)
CC_NSOvsNSI=read.table("~/Documents/JAMES Broken Computer/Guam/Raymundo/Reanalysis/In_Situ/MWU_CC_NS_OvsNS_I_GO.csv", header=T)
NSOvsNSI=rbind(BP_NSOvsNSI,MF_NSOvsNSI,CC_NSOvsNSI)

BP_FSOvsFSI=read.table("~/Documents/JAMES Broken Computer/Guam/Raymundo/Reanalysis/In_Situ/MWU_BP_FSOvsFSI_GO.csv", header=T)
MF_FSOvsFSI=read.table("~/Documents/JAMES Broken Computer/Guam/Raymundo/Reanalysis/In_Situ/MWU_MF_FSOvsFSI_GO.csv", header=T)
CC_FSOvsFSI=read.table("~/Documents/JAMES Broken Computer/Guam/Raymundo/Reanalysis/In_Situ/MWU_CC_FSOvsFSI_GO.csv", header=T)
FSOvsFSI=rbind(BP_FSOvsFSI,MF_FSOvsFSI,CC_FSOvsFSI)

BP_FSvsNS=read.table("~/Documents/JAMES Broken Computer/Guam/Raymundo/Reanalysis/In_Situ/MWU_BP_FSvsNS_GO.csv",header=T)
MF_FSvsNS=read.table("~/Documents/JAMES Broken Computer/Guam/Raymundo/Reanalysis/In_Situ/MWU_MF_FSvsNS_GO.csv",header=T)
CC_FSvsNS=read.table("~/Documents/JAMES Broken Computer/Guam/Raymundo/Reanalysis/In_Situ/MWU_CC_FSvsNS_GO.csv",header=T)
#FSvsNS=rbind(BP_FSvsNS,MF_FSvsNS,CC_FSvsNS)
#renamed
insitu_HFvsLF=rbind(BP_FSvsNS,MF_FSvsNS,CC_FSvsNS)


BP_FSvsNS_B=read.table("~/Documents/JAMES Broken Computer/Guam/Raymundo/Reanalysis/In_Situ/MWU_BP_FSvsNS_B_GO.csv",header=T)
MF_FSvsNS_B=read.table("~/Documents/JAMES Broken Computer/Guam/Raymundo/Reanalysis/In_Situ/MWU_MF_FSvsNS_B_GO.csv",header=T)
CC_FSvsNS_B=read.table("~/Documents/JAMES Broken Computer/Guam/Raymundo/Reanalysis/In_Situ/MWU_CC_FSvsNS_B_GO.csv",header=T)
#FSvsNS_B=rbind(BP_FSvsNS_B,MF_FSvsNS_B,CC_FSvsNS_B)
#renamed
insitu_HFvsLF_B=rbind(BP_FSvsNS_B,MF_FSvsNS_B,CC_FSvsNS_B)

BP_FSBvsFSU=read.table("~/Documents/JAMES Broken Computer/Guam/Raymundo/Reanalysis/In_Situ/MWU_BP_BvsU_FS_GO.csv",header=T)
MF_FSBvsFSU=read.table("~/Documents/JAMES Broken Computer/Guam/Raymundo/Reanalysis/In_Situ/MWU_MF_BvsU_FS_GO.csv",header=T)
CC_FSBvsFSU=read.table("~/Documents/JAMES Broken Computer/Guam/Raymundo/Reanalysis/In_Situ/MWU_CC_BvsU_FS_GO.csv",header=T)
#FSBvsFSU=rbind(BP_FSBvsFSU,MF_FSBvsFSU,CC_FSBvsFSU)
#renamed
insitu_HFBvsHFU=rbind(BP_FSBvsFSU,MF_FSBvsFSU,CC_FSBvsFSU)

BP_NSBvsNSU=read.table("~/Documents/JAMES Broken Computer/Guam/Raymundo/Reanalysis/In_Situ/MWU_BP_BvsU_NS_GO.csv",header=T)
MF_NSBvsNSU=read.table("~/Documents/JAMES Broken Computer/Guam/Raymundo/Reanalysis/In_Situ/MWU_MF_BvsU_NS_GO.csv",header=T)
CC_NSBvsNSU=read.table("~/Documents/JAMES Broken Computer/Guam/Raymundo/Reanalysis/In_Situ/MWU_CC_BvsU_NS_GO.csv",header=T)
#NSBvsNSU=rbind(BP_NSBvsNSU,MF_NSBvsNSU,CC_NSBvsNSU)
#renamed
insitu_LFBvsLFU=rbind(BP_NSBvsNSU,MF_NSBvsNSU,CC_NSBvsNSU)

library(dplyr)
library(ggplot2)
library(ggrepel)


##CORRELATIONS##

##rplot plots all GO terms
##GGplot plots all GO terms plus add labels for GO terms sig in BOTH comparisons

#PLOT only looking at sig GOs under either of the two treatments (see below for all GOs)
test2=function(a,b) {
  goods=intersect(a$term,b$term);print(length(goods)); 
  data1=a[a$term %in% goods,];
  data2=b[b$term %in% goods,]; 
  ress=merge(data1,data2,by="term");
  # GO terms signifcant in both datasets
  sigs=(ress$p.adj.x<=0.1 & ress$p.adj.y<=0.1);
  eithersigs=(ress$p.adj.x<=0.1 | ress$p.adj.y<=0.1)
  EITHERSIGS=ress[eithersigs,]
  SIGS=ress[sigs,]
  #shapiro.test(EITHERSIGS$delta.rank.x)
  #shapiro.test(EITHERSIGS$delta.rank.y)
  COR=cor.test(EITHERSIGS$delta.rank.x,EITHERSIGS$delta.rank.y,method="pearson")
  #ress.lm = lm(delta.rank.x ~ delta.rank.y, data=EITHERSIGS);
  name1 <- deparse(substitute(a));
  name2<-deparse(substitute(b));
  plot(delta.rank.x~delta.rank.y,EITHERSIGS,
       main=paste(name1,"vs",name2),
       sub=paste("pvalue=",COR[["p.value"]]),
       xlab="GBM", ylab="GE",mgp=c(2.3,1,0))
  abline(v=0,lty=3)
  abline(h=0,lty=3)
  abline(lm(ress$delta.rank.x ~ ress$delta.rank.y))
  ress.lm = lm(delta.rank.x ~ delta.rank.y, data=EITHERSIGS)
  print(COR);
  pval=signif(COR[["p.value"]],3);
  r=signif(COR[["estimate"]],2);
  ggplot(EITHERSIGS, aes(delta.rank.y, delta.rank.x))+
    geom_point()+geom_hline(yintercept = 0)+
    geom_vline(xintercept=0)+
    stat_smooth(method = "lm", col = "red")+
    xlab(paste(name2))+ylab(paste(name1))+ggtitle(paste(name1,"vs",name2), subtitle=paste("p=",pval,  "r=",r))+
    # geom_label_repel(aes(label = SIGS$name.x),
    #                  box.padding   = 0.15, 
    #                  point.padding = 1,
    #                  segment.color = 'blue',data=SIGS, segment.size = 1.0) +
    theme_classic()+
    theme(text = element_text(size = 20))
}





#Comparing ambient HF/LF comparison to in situ heat stress(aka before bleaching) HF/LF
test(HFNHvsLFNH, insitu_HFvsLF)
#Comparing heat stress HF/LF comparison to in situ bleaching HF/LF
test(HFHvsLFH,insitu_HFvsLF_B)
#Comparing heat stress HF/LF to in situ heat stress
test(HFHvsLFH,insitu_HFvsLF)

#2
#Comparing ambient HF/LF comparison to in situ heat stress(aka before bleaching) HF/LF
test2(HFNHvsLFNH, insitu_HFvsLF)
#Comparing heat stress HF/LF comparison to in situ bleaching HF/LF
# test2(HFHvsLFH,insitu_HFvsLF_B)
# #Comparing heat stress HF/LF to in situ heat stress
# test2(HFHvsLFH,insitu_HFvsLF)
# ## Looking at within site OvsI 
# test2(HFHvsLFH, FSOvsFSI)
# test2(HFHvsLFH,NSOvsNSI)
# test2(HFNHvsLFNH, FSOvsFSI)
# test2(HFNHvsLFNH, NSOvsNSI)
# test2(FSOvsFSI,NSOvsNSI)
#  
#looking at FS bleached vs unbleached overlap with
#Looking at NS bleached vs unbleached
test2(insitu_HFBvsHFU,insitu_LFBvsLFU)





####EXTRA#########################################################
##REGRESSIONS##

test=function(a,b) {
  goods=intersect(a$term,b$term);print(length(goods)); 
  data1=a[a$term %in% goods,];
  data2=b[b$term %in% goods,]; 
  ress=merge(data1,data2,by="term");
  # GO terms signifcant in both datasets
  sigs=(ress$p.adj.x<=0.1 & ress$p.adj.y<=0.1);
  SIGS=ress[sigs,]
  ress.lm = lm(delta.rank.x ~ delta.rank.y, data=ress);
  name1 <- deparse(substitute(a));
  name2<-deparse(substitute(b));
  plot(delta.rank.x~delta.rank.y,ress,
       main=paste(name1,"vs",name2),
       sub=paste("pvalue=",summary(ress.lm)$coefficients[2,4]),
       xlab="GBM", ylab="GE",mgp=c(2.3,1,0))
  abline(v=0,lty=3)
  abline(h=0,lty=3)
  abline(lm(ress$delta.rank.x ~ ress$delta.rank.y))
  ress.lm = lm(delta.rank.x ~ delta.rank.y, data=ress)
  print(summary(ress.lm));
  ggplot(ress, aes(delta.rank.y, delta.rank.x))+
    geom_point()+geom_hline(yintercept = 0)+
    geom_vline(xintercept=0)+
    stat_smooth(method = "lm", col = "red")+
    xlab("GBM")+ylab("GE")+ggtitle(paste(name1,"vs",name2), subtitle=paste("pvalue=",summary(ress.lm)$coefficients[2,4]))+
    geom_label_repel(aes(label = SIGS$name.x),
                     box.padding   = 0.15, 
                     point.padding = 1,
                     segment.color = 'blue',data=SIGS, segment.size = 1.0) +
    theme_classic()
  
}

##rplot plots all GO terms
##GGplot plots all GO terms plus add labels for GO terms sig in BOTH comparisons

#PLOT only looking at sig GOs under either of the two treatments (previous is all GOs)
test2=function(a,b) {
  goods=intersect(a$term,b$term);print(length(goods)); 
  data1=a[a$term %in% goods,];
  data2=b[b$term %in% goods,]; 
  ress=merge(data1,data2,by="term");
  # GO terms signifcant in both datasets
  sigs=(ress$p.adj.x<=0.1 & ress$p.adj.y<=0.1);
  eithersigs=(ress$p.adj.x<=0.1 | ress$p.adj.y<=0.1)
  EITHERSIGS=ress[eithersigs,]
  SIGS=ress[sigs,]
  ress.lm = lm(delta.rank.x ~ delta.rank.y, data=EITHERSIGS);
  name1 <- deparse(substitute(a));
  name2<-deparse(substitute(b));
  plot(delta.rank.x~delta.rank.y,EITHERSIGS,
       main=paste(name1,"vs",name2),
       sub=paste("pvalue=",summary(ress.lm)$coefficients[2,4]),
       xlab="GBM", ylab="GE",mgp=c(2.3,1,0))
  abline(v=0,lty=3)
  abline(h=0,lty=3)
  abline(lm(ress$delta.rank.x ~ ress$delta.rank.y))
  ress.lm = lm(delta.rank.x ~ delta.rank.y, data=EITHERSIGS)
  print(summary(ress.lm));
  pval=signif(summary(ress.lm)$coefficients[2,4],3);
  r2=signif(summary(ress.lm)$r.squared,2);
  ggplot(EITHERSIGS, aes(delta.rank.y, delta.rank.x))+
    geom_point()+geom_hline(yintercept = 0)+
    geom_vline(xintercept=0)+
    stat_smooth(method = "lm", col = "red")+
    xlab(paste(name2))+ylab(paste(name1))+ggtitle(paste(name1,"vs",name2), subtitle=paste("p=",pval,  "R²=",r2))+
    # geom_label_repel(aes(label = SIGS$name.x),
    #                  box.padding   = 0.15, 
    #                  point.padding = 1,
    #                  segment.color = 'blue',data=SIGS, segment.size = 1.0) +
    theme_classic()
}

###Same code from above, but broken up ##


#############MF############# MF###################
setwd("~/Documents/JAMES Broken Computer/Guam/Raymundo/Reanalysis/Ex_Situ")
data1=read.table("MWU_MF_BvsD_GO.csv",header=T) #ambient temp comp
data1=read.table("MWU_MF_AvsC_GO.csv",header=T) #hs comp
#data1=read.table("MWU_BP_CvsD_GO.csv", header=T) #low flow comp
#data1=read.table("MWU_BP_AvsB_GO.csv", header=T) #high flow comp
setwd("~/Documents/JAMES Broken Computer/Guam/Raymundo/Reanalysis/In_Situ/")
data2=read.table("MWU_MF_FSvsNS_GO.csv",header=T)
data2=read.table("MWU_BP_FSvsNS_B_GO.csv",header=T)
#data2=read.table("MWU_BP_BvsU_NS_GO.csv") # had nothing in common with ex situ
goods=intersect(data1$term,data2$term)
#goods=unique(as.character(c(data1$term[data1$p.adj<=0.1],data2$term[data2$p.adj<=0.1])))
length(goods)

data1=data1[data1$term %in% goods,]
data2=data2[data2$term %in% goods,]

# all overlapping GO terms
ress=merge(data1,data2,by="term")
plot(delta.rank.x~delta.rank.y,ress,xlab="GBM", ylab="GE",mgp=c(2.3,1,0))
abline(v=0,lty=3)
abline(h=0,lty=3)

# GO terms highly signifcant in any of the two datasets
sigs=(ress$p.adj.x<=0.01 | ress$p.adj.y<=0.01)
sum(sigs) # 71
plot(delta.rank.x~delta.rank.y,ress[sigs,],xlab="GBM", ylab="GE",mgp=c(2.3,1,0))
abline(v=0,lty=3)
abline(h=0,lty=3)

# GO terms signifcant in both datasets
sigs=(ress$p.adj.x<=0.1 & ress$p.adj.y<=0.1)
sum(sigs) # 20
plot(delta.rank.x~delta.rank.y,ress[sigs,],xlab="GBM", ylab="GE",mgp=c(2.3,1,0))
abline(v=0,lty=3)
abline(h=0,lty=3)

MFpre=ress[sigs,]
# MFpost had none in common
########################## CC###################
setwd("~/Documents/JAMES Broken Computer/Guam/Raymundo/Reanalysis/Ex_Situ")
data1=read.table("MWU_CC_BvsD_GO.csv",header=T) #ambient temp comp
data1=read.table("MWU_CC_AvsC_GO.csv",header=T) #hs comp
#data1=read.table("MWU_BP_CvsD_GO.csv", header=T) #low flow comp
#data1=read.table("MWU_BP_AvsB_GO.csv", header=T) #high flow comp
setwd("~/Documents/JAMES Broken Computer/Guam/Raymundo/Reanalysis/In_Situ/")
data2=read.table("MWU_CC_FSvsNS_GO.csv",header=T)
data2=read.table("MWU_CC_FSvsNS_B_GO.csv",header=T)
#data2=read.table("MWU_BP_BvsU_NS_GO.csv") # had nothing in common with ex situ
goods=intersect(data1$term,data2$term)
#goods=unique(as.character(c(data1$term[data1$p.adj<=0.1],data2$term[data2$p.adj<=0.1])))
length(goods)

data1=data1[data1$term %in% goods,]
data2=data2[data2$term %in% goods,]

# all overlapping GO terms
ress=merge(data1,data2,by="term")
plot(delta.rank.x~delta.rank.y,ress,xlab="GBM", ylab="GE",mgp=c(2.3,1,0))
abline(v=0,lty=3)
abline(h=0,lty=3)

# GO terms highly signifcant in any of the two datasets
sigs=(ress$p.adj.x<=0.01 | ress$p.adj.y<=0.01)
sum(sigs) # 71
plot(delta.rank.x~delta.rank.y,ress[sigs,],xlab="GBM", ylab="GE",mgp=c(2.3,1,0))
abline(v=0,lty=3)
abline(h=0,lty=3)

# GO terms signifcant in both datasets
sigs=(ress$p.adj.x<=0.1 & ress$p.adj.y<=0.1)
sum(sigs) # 20
plot(delta.rank.x~delta.rank.y,ress[sigs,],xlab="GBM", ylab="GE",mgp=c(2.3,1,0))
abline(v=0,lty=3)
abline(h=0,lty=3)

CCpre=ress[sigs,]
CCpost=ress[sigs,] #none


#########Graphing
AllBothSigsPre= rbind(BPpre, MFpre, CCpre)
library(ggplot2)
library(ggrepel)
f <- ggplot(AllBothSigsPre, aes(delta.rank.x, delta.rank.y))+
  geom_point()+geom_hline(yintercept = 0)+
  geom_vline(xintercept=0)+
  xlab("GBM")+ylab("GE")+ggtitle("Delta rank 'pre-stress'")

AllBothSigsPrePlot=f +
  geom_label_repel(aes(label = name.x),
                   box.padding   = 0.15, 
                   point.padding = 0.25,
                   segment.color = 'grey50') +
  theme_classic()
AllBothSigsPrePlot
dev.copy(png,'AllBothSigsPrePlot.png')

dev.off()
#######Post
AllBothSigsPost= rbind(BPpost)
library(ggplot2)
library(ggrepel)
f <- ggplot(AllBothSigsPost, aes(delta.rank.x, delta.rank.y))+
  geom_point()+geom_hline(yintercept = 0)+
  geom_vline(xintercept=0)+
  xlab("GBM")+ylab("GE")+ggtitle("Delta rank 'stressed'")

AllBothSigsPostPlot=f +
  geom_label_repel(aes(label = name.x),
                   box.padding   = 0.15, 
                   point.padding = 0.25,
                   segment.color = 'grey50') +
  theme_classic()
AllBothSigsPostPlot
dev.copy(png,'AllBothSigsPrePlot.png')

dev.off()
