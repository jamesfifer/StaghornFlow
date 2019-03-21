#!/usr/bin/env python
from Bio import SeqIO

fasta_file = "/media/RAID/james/Staghorn_GeneExpression/Ref_Transcriptome_All/Filtered_Transcriptome/Trinity-GG.final.fasta" # Input fasta file
wanted_file = "/media/RAID/james/Staghorn_GeneExpression/Ref_Transcriptome_All/Annotation/BlastxAnnotated.ID.txt" # Input interesting sequence IDs, one per line
result_file = "/media/RAID/james/Staghorn_GeneExpression/Ref_Transcriptome_All/Annotation/Trinity-GG.annotated.fasta" # Output fasta file

wanted = set()
with open(wanted_file) as f:
    for line in f:
        line = line.strip()
        if line != "":
            wanted.add(line)

fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
with open(result_file, "w") as f:
    for seq in fasta_sequences:
        if seq.id in wanted:
            SeqIO.write([seq], f, "fasta")
