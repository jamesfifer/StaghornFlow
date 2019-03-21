#Take DE subsets with IDs from relevant comparisons

#awk '{print $1}' ../DESeq*/RSEM.gene.counts.matrix.Cond_B_vs_Cond_D.DESeq2.DE_results.P0.05_C1.DE.subset > ./HFvLF.DESeq2.txt

#awk '{print $1}' ../edgeR*/RSEM.gene.counts.matrix.Cond_B_vs_Cond_D.edgeR.DE_results.P0.05_C1.DE.subset > ./HFvLF.edgeR.txt

# Find intersection

#cat ./HFvLF.DESeq2.txt | grep -F -f HFvLF.edgeR.txt >HFvLF.DE.genes.txt
#grep -F -f HFvLF.DE.genes.txt ../DESeq*/RSEM.gene.counts.matrix.Cond_B_vs_Cond_D.DESeq2.DE_results.P0.05_C1.Cond_B-UP.subset | awk '{print $1}'  > HF.UP.DE.genes.txt
#grep -F -f HFvLF.DE.genes.txt ../DESeq*/RSEM.gene.counts.matrix.Cond_B_vs_Cond_D.DESeq2.DE_results.P0.05_C1.Cond_D-UP.subset | awk '{print $1}'  > HF.DOWN.DE.genes.txt

#Obtain annotation 
#grep -F -f HFvLF.DE.genes.txt /media/RAID/james/Staghorn_GeneExpression/Ref_Transcriptome_All/Annotation/BlastxAnnotated.eValueCutoff.txt >HFvLF.annotated.txt
#Look for candidate genes
#grep -F -f HFvLF.DE.genes.txt /media/RAID/james/Staghorn_GeneExpression/Ref_Transcriptome_All/Candidate/Candidate.genes.txt >HFvLF.candidate.txt
#See what candidate genes are upregulated/downregulated

#See what other annotated genes are upregulated/downregulated

#
#awk '{print $1}' HFvLF.candidate.txt > ID; sed 's/...$//' ID > newID; \
#grep -F -f newID ../DESeq2*/RSEM.gene.counts.matrix.Cond_B_vs_Cond_D.DESeq2.DE_results.P0.05_C1.Cond_B-UP.subset\
# | awk '{print $1,$7}' > Up; grep -F -f ID ../DESeq2*/RSEM.gene.counts.matrix.Cond_B_vs_Cond_D.DESeq2.DE_results.P0.05_C1.Cond_D-UP.subset\
# | awk '{print $1,$7}' > Down; cat Up Down > HFvLF.candidate.fc.txt; rm Up Down ID newID;

#grep -F -f HFvLF.annotated.ID.txt ../DESeq2*/RSEM.gene.counts.matrix.Cond_B_vs_Cond_D.DESeq2.DE_results.P0.05_C1.Cond_B-UP.subset\
# | awk '{print $1,$7}' > Up; grep -F -f HFvLF.annotated.ID.txt ../DESeq2*/RSEM.gene.counts.matrix.Cond_B_vs_Cond_D.DESeq2.DE_results.P0.05_C1.Cond_D-UP.subset\
# | awk '{print $1,$7}' > Down; cat Up Down > HFvLF.annotation.fc.txt; rm Up Down ID newID;

grep -F -f HF.DOWN.DE.genes.txt HFvLF.annotated.txt > HFvLF.Down.annotated.txt
grep -F -f HF.UP.DE.genes.txt HFvLF.annotated.txt > HFvLF.Up.annotated.txt 
