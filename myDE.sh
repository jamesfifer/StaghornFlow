#!/bin/bash
#below code most be run while in directory with samples  otherwise you will get unable to locate sample__
#Trinity recommends running align and estimate abundance first without sample file to prep reference and then with sample file without -prep reference command. (Prep the reference and then run on samples only). 
#/usr/local/src/trinityrnaseq/util/align_and_estimate_abundance.pl --transcripts 
/media/RAID/james/Staghorn_GeneExpression/Ref_Transcriptome_All/Annotation/Trinity-GG.final.fasta  --est_method RSEM --aln_method bowtie2 
--prep_reference --trinity_mode 
#/usr/local/src/trinityrnaseq/util/align_and_estimate_abundance.pl --transcripts /media/RAID/james/Staghorn_GeneExpression/Ref_Transcriptome_All/Filtered_Transcriptome/Trinity-GG.final.fasta --est_method RSEM --aln_method bowtie2 --output_dir /media/RAID/james/Staghorn_GeneExpression/Ex_Situ/DE_Analysis_All/ --trinity_mode --samples /media/RAID/james/Staghorn_GeneExpression/Ex_Situ/DE_Analysis_All/new.samples.txt --seqType fq

#/usr/local/src/trinityrnaseq/util/abundance_estimates_to_matrix.pl \
#	--est_method RSEM \
#	--gene_trans_map /media/RAID/james/Staghorn_GeneExpression/Ref_Transcriptome_All/Filtered_Transcriptome/Trinity-GG.final.fasta.gene_trans_map \
#	--name_sample_by_basedir ./Cond_A_rep1/RSEM.isoforms.results ./Cond_A_rep2/RSEM.isoforms.results ./Cond_A_rep3/RSEM.isoforms.results ./Cond_A_rep4/RSEM.isoforms.results \
#	./Cond_B_rep1/RSEM.isoforms.results ./Cond_B_rep2/RSEM.isoforms.results ./Cond_B_rep3/RSEM.isoforms.results ./Cond_B_rep4/RSEM.isoforms.results \
#	./Cond_C_rep1/RSEM.isoforms.results ./Cond_C_rep2/RSEM.isoforms.results ./Cond_C_rep3/RSEM.isoforms.results ./Cond_C_rep4/RSEM.isoforms.results \
#	./Cond_D_rep1/RSEM.isoforms.results ./Cond_D_rep3/RSEM.isoforms.results ./Cond_D_rep4/RSEM.isoforms.results \


 # R
#  source("http://bioconductor.org/biocLite.R")
#  biocLite('edgeR')
#  biocLite('DESeq2')
#  biocLite('ctc')
#  biocLite('Biobase')
#  install.packages('gplots')
#  install.packages('ape')
#--matrix|m RSEM.gene.counts.matrix \


#QUALITY CHECK (Trinity recommends doing this before DE to identify potential outliers/batch effects, but don't read too much into it) Examine replicates 
#/usr/local/src/trinityrnaseq/Analysis/DifferentialExpression/PtR \
# -m RSEM.gene.counts.matrix -s samples.txt \
# --log2 --compare_replicates 
#Examine all samples:
#/usr/local/src/trinityrnaseq/Analysis/DifferentialExpression/PtR \
#      -m RSEM.gene.counts.matrix -s samples.txt \
#      --log2 --sample_cor_matrix

#Could also run the below code on isoform files to get transcript DE analysis. 
#/usr/local/src/trinityrnaseq/util/..//Analysis/DifferentialExpression/run_DE_analysis.pl  \
#--method DESeq2 \
#--matrix /media/RAID/james/Staghorn_GeneExpression/Ex_Situ/DE_Analysis_All/RSEM.gene.counts.matrix \
#--samples_file /media/RAID/james/Staghorn_GeneExpression/Ex_Situ/DE_Analysis_All/new.samples.txt 

#/usr/local/src/trinityrnaseq/util/..//Analysis/DifferentialExpression/run_DE_analysis.pl  \
#--method edgeR \
#--matrix /media/RAID/james/Staghorn_GeneExpression/Ex_Situ/DE_Analysis_All/RSEM.gene.counts.matrix \
#--samples_file /media/RAID/james/Staghorn_GeneExpression/Ex_Situ/DE_Analysis_All/new.samples.txt \

#/usr/local/src/trinityrnaseq/Analysis/DifferentialExpression/PtR \
#       -m RSEM.gene.counts.matrix -s samples.txt\
#       --log2 --CPM --prin_comp 3

#
#/usr/local/src/trinityrnaseq/Analysis/DifferentialExpression/analyze_diff_expr.pl \
#      --matrix ../RSEM.gene.TMM.EXPR.matrix \
#      --samples ../samples.txt \
#      -P 1e-3 -C 2

#C refers to fold change (2^C fold change)
#Need to run in DESeq or edgeR directory
#Run with p value .001, .01 and .05

#/usr/local/src/trinityrnaseq/util/..//Analysis/DifferentialExpression/analyze_diff_expr.pl \
#--matrix ../RSEM.gene.TMM.EXPR.matrix --P 0.05 --C 1 --samples ../new.samples.txt

#For full pipeline example goto 
#/usr/local/src/trinityrnaseq/trinity_ext_sample_data/test_full_edgeR_pipeline

