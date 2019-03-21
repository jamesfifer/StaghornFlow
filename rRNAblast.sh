
#make database

#makeblastdb -in ./arb-silva.de_2018-02-07_id497276_tax_silva.fasta -dbtype nucl

#nice -6 blastn -query ./Trinity-GG.NCBI.fasta -db /media/RAID/james/Genomes/rRNA/arb-silva.de_2018-02-07_id497276_tax_silva.fasta -out rRNA_Filt.blast.txt -evalue 1e-8 -outfmt 6 -max_target_seqs 1 -num_threads 40 

#Scripts to create fasta file from blast hit file:

#awk '{print $1}' rRNA_Coral.blast.txt  > rRNAHits.txt
#cat rRNA_Coral.blast.txt  | grep -f rRNAHits.txt > rRNAHits.fasta

#script to remove rRNA from GG Transcriptome:
#id=$(cat rRNA_Filt.ID.txt) ;
#for gene in $id; 
#do 
#awk '/'$gene'/{flag=1;print $0;next}/^>/{flag=0}flag' Trinity-GG.Filt.fasta >> Trinity-GG.FINAL.fasta ;
#done

# use below after obtaining file with gene names only
#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from Bio import SeqIO

fasta_file = sys.argv[1]  # Input fasta file
remove_file = sys.argv[2] # Input file to remove, one gene name per line
result_file = sys.argv[3] # Output fasta file

remove = set()
with open(remove_file) as f:
    for line in f:
        line = line.strip()
        if line != "":
            remove.add(line)

fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')

with open(result_file, "w") as f:
    for seq in fasta_sequences:
        nuc = seq.seq.tostring()
        if nuc not in remove and len(nuc) > 0:
            SeqIO.write([seq], f, "fasta")
