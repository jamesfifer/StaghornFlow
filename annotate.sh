#!/bin/bash

#Creates the necessary .pep file:
#~/Transcriptome_Tools/TransDecoder-TransDecoder-v5.0.2/TransDecoder.LongOrfs -t /media/RAID/james/Staghorn_GeneExpression/Ref_Transcriptome_All/Filtered_Transcriptome/Trinity-GG.final.fasta
#nice -n 6 ~/Transcriptome_Tools/TransDecoder-TransDecoder-v5.0.2/TransDecoder.Predict -t 
/media/RAID/james/Staghorn_GeneExpression/Ref_Transcriptome_All/Trinity-GG.final.fasta

#Download the uniprot database you want with GO terms.
# can do this with Kegg as well
makeblastdb -in uniprot-cnidariaGO.fasta -dbtype prot

#Blast
nice -n 6 blastx -query /media/RAID/james/Staghorn_GeneExpression/Ref_Transcriptome_All/Filtered_Transcriptome/Trinity-GG.final.fasta -db ./uniprot-cnidariakegg.fasta -max_hsps 1 -num_threads 40 -max_target_seqs 1 -outfmt 6 > keggblastx.outfmt6
#We only used blastx for our annotations, but you could theoretically do blastp as well. 


# the annotation output from trinotate is sparse and unhelpful, use the script ParseAnnotations.py instead once you have the blast 
outputs





