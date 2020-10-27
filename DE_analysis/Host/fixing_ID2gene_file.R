#fixing ID2gene file
library(dplyr)
library(tidyr)
annotation<-read.csv("ID2gene.noGO.txt", sep='\t', header=F)
new<-annotation %>%
  separate(V1, c("ID", "gene"), "_i[0-9]")
new1<-new[!duplicated(new$ID),]

write.table(new1, file="ID2gene.tab", quote=F, sep="\t")
