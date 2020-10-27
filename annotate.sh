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

####Alternatively use plast (we found no real differences between plast and blast searches, but plast is exponentially faster (like 100x faster)

#/projectnb/incrna/Plast/build/bin/plast -p plastx -i Sym.fasta -d ./uniprot_sprot.fasta -max-hit-per-query 1 -max-hsp-per-hit 1 -e 1e-4 -o ./uniprot.plast.out

# make a file that only has the unique results (save lowest e-value for each transcript)
#cat ./uniprot.plast.out | sort -k1,1 -k11,11g | awk '!seen[$1]++' > uniq_plastresults.out

#grep ">" Sym.fasta| perl -pe 's/>comp(\d+)(\S+)+/comp$1$2\tisogroup$1/' > Sym_seq2iso.tab

#grep ">" Sym.fasta | perl -pe 's/>comp(\d+)(\S+)\s.+/comp$1$2\tisogroup$1/' > sym_seq2iso.tab
#cat Sym.fasta | perl -pe 's/>comp(\d+)(\S+)+/>comp$1$2 gene=isogroup$1/' > sym_iso.fasta
#cat uniq_plastresults.out | sed 's/comp/isogroup/' |  sed 's/_c[0-9]_seq[0-9]//'> uniq_plastresults_isogroup.out
#grep ">" uniprot_sprot.fasta | awk '{print $1}' | sed 's/>sp[|]//' | sed 's/[|].*//' > geneIDs
#grep ">" uniprot_sprot.fasta | awk '{$1="";print $0}' > geneDescriptions
#paste geneIDs geneDescriptions > longGeneNames
#cat uniq_plastresults_isogroup.out  | awk {'print $1"\t"$2'} > iso2gene.tab
#cat iso2gene.tab | awk -F '[|]' {'print $1"\t"$2'} | cut -f1,3 > iso2shortgene.tab
#join -1 2 -2 1 -e "noMatch" -t $'\t' <(sort -k2 iso2shortgene.tab) <(sort -k1 ./longGeneNames) | awk 'BEGIN {FS="\t"; OFS="\t"} {print $2, $1, $3}' > Sym_iso2geneDescription.tab
#wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz
#12/18/19 version 16.7GB
#gunzip *.gz


#join -1 2 -2 1 -e "noMatch" -t $'\t' <(sort -k2 iso2shortgene.tab) <(sort -k1 ./idmapping_selected.tab) | awk 'BEGIN {FS="\t"; OFS="\t"} {print $2, $1, $8}' > Sym_iso2go.tab

awk '{print $1}' sym_iso2go.tab | cut -f1,2,3,4,5 -d'_' | paste - sym_iso2go.tab | cut -f1,4- > Sym_ID2GO.tab
awk '$2 !~ "noMatch"' Sym_ID2GO.tab > Sym_ID2GO.tab.tmp
#then use nrify_GOtable.pl
./nrify_GOtable.pl Sym_ID2GO.tab.tmp > Sym_ID2GO.tab
sed 's/ //1g' Sym_ID2GO.tab > Sym_ID2GO.tab2


awk '{print $1}' Sym_iso2geneDescription.tab | cut -f1,2,3,4,5 -d'_' | paste - Sym_iso2geneDescription.tab | cut -f1,4- > Sym_ID2Description.tab
#some isogroups had different descriptions so only kept first description when there were duplicates
#if you wanted to double check if this would be deleting an interesting description you can first check what they look like
awk '
  !count[$1]++ {save[$1] = $0; next}
  count[$1] == 2 {
    print save[$1]
    delete save[$1]
  }
  {print}' Sym_ID2Description.tab
#otherwise remove
sort -u -k1,1 Sym_ID2Description.tab >Sym_ID2Description.tab.tmp
#
mv Sym_ID2Description.tab.tmp Sym_ID2Description.tab








