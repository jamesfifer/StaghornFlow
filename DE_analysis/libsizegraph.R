setwd("~/Documents/JAMES Broken Computer/Guam/Raymundo/Reanalysis")
exsampdata=read.table("~/Documents/JAMES Broken Computer/Guam/Raymundo/Reanalysis/Ex_Situ/samples.txt", header=F)
treat=c( "HFH", "HFH", "HFH" , "HFH", "HFNH","HFNH", "HFNH", "HFNH", "LFH","LFH","LFH","LFH","LFNH","LFNH","LFNH")
colony=as.factor (c( "b", "a", "c", "d", "a","b", "c", "d", "b","a","c","d","b","a","c"))
exp=as.factor(c("ex", "ex", "ex" , "ex", "ex","ex", "ex", "ex", "ex","ex","ex","ex","ex","ex","ex"))
g=data.frame(treat, colony,exp)
g$cond=paste(g$treat,g$colony, sep=".")

ex=cbind(exsampdata, g)[c(3,7,8)]

insampdata=read.table("~/Documents/JAMES Broken Computer/Guam/Raymundo/Reanalysis/In_Situ/samples.txt", header=F)
treat=c( "HF", "HF", "HF" , "LF", "LF","LF", "HF", "HF", "HF","HF","HF","HF","LF","LF","LF","LF","LF","LF")
colony=as.factor (c( "Col_1", "Col_2", "Col_3", "Col_4", "Col_5","Col_6", "Col_1", "Col_2", "Col_3","Col_1","Col_2","Col_3","Col_4","Col_5","Col_6","Col_4","Col_5","Col_6"))
O.I=c("O", "O", "O","O", "O", "O","O", "O", "O","I","I","I","O","O","O","I","I","I")
bleach.state= c( "B", "B", "B" , "B", "B","B", "HS", "HS", "HS","HS","HS","HS","HS","HS","HS","HS","HS","HS")
exp=c( "in", "in", "in" , "in", "in","in", "in", "in", "in","in","in","in","in","in","in","in","in","in")

g=data.frame(treat,colony,O.I, bleach.state,exp)
g$cond=paste(g$treat,g$O.I,g$bleach.state,g$colony, sep=".")

IN=cbind(insampdata, g)[c(3,9,10)]


all=rbind(ex,IN)
colnames(all)[1] = "V2"


lib=read.table("~/Documents/JAMES Broken Computer/Guam/Raymundo/Reanalysis/librarysizes.txt", sep='\t', header=F)
newlib=lib[!duplicated(lib[,c('V1')]),]
newlib$V3<-gsub("_R1_Combined.fq", "",newlib$V2)
newlib=newlib[-c(20,35),,drop=FALSE] 

newcd<-transform(merge(newlib,all,by="V2"))

insitu=filter(newcd, exp=="in")
insitu <- insitu[order(insitu$cond),]
exsitu=filter(newcd,exp=="ex")
exsitu <- exsitu[order(exsitu$cond),]


par(mar=c(5,6,4,1)+2)
barplot(insitu$V1, main="Library size paired end reads in-situ",
        col=c("grey"),
        names.arg = insitu$cond,las=2, ylab="lib size (in millions)")

barplot(exsitu$V1, main="Library size paired end reads ex-situ",
        col=c("grey"),
        names.arg = exsitu$cond,las=2, ylab="lib size (in millions)")

################
mappingrateH=read.table(file = "~/Documents/JAMES Broken Computer/Guam/Raymundo/Reanalysis/mapping.host.sum.stats.txt")
mappingrateH$V3<-gsub("/align_summary.txt", "",mappingrateH$V1)
mappingrateH$V3<-gsub("./", "",mappingrateH$V3)
#fix bfs names
mappingrateH[3, 3] <- "BF1_S12"
mappingrateH[12, 3] <- "BF9_S14"
mappingrateH[22,3]<- "BF5_S12"
mappingrateH[17,3]="BN5_S10"
mappingrateH$type="Host"
mappingrateH <- mappingrateH[order(mappingrateH$V3),]
mappingrateH<-transform(merge(mappingrateH,newcd,by="V3"))



mappingrateS=read.table(file = "~/Documents/JAMES Broken Computer/Guam/Raymundo/Reanalysis/mapping.sym.sum.stats.txt")
mappingrateS$V3<-gsub("/align_summary.txt", "",mappingrateS$V1)
mappingrateS$V3<-gsub("./", "",mappingrateS$V3)
mappingrateS$type="Sym"
mappingrateS <- mappingrateS[order(mappingrateS$V3),]
mappingrateS<-transform(merge(mappingrateS,newcd,by="V3"))

maprates=rbind(mappingrateH,mappingrateS)
maprates$V2.x=as.numeric(sub("%", "", maprates$V2.x))

library(ggplot2)
p <-ggplot(maprates, aes(cond, V2.x))
p +geom_bar(stat = "identity", aes(fill = type)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+
  ylab("% mapped")+xlab(NULL)+ggtitle("Mapping rates")




