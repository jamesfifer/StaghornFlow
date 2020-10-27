#!/bin/bash
#below code most be run while in directory with samples  otherwise you will get unable to locate sample__
#Trinity recommends running align and estimate abundance first without sample file to prep reference and then with sample file without -prep reference command. (Prep the reference and then run on samples only).
module load bowtie2/2.3.4.1
module load rsem/1.3.1
module load htslib/1.9
module load samtools/1.9
#prep reference first
#/projectnb/davieslab/jfifer/Sym_Flow/trinityrnaseq-v2.9.0/util/align_and_estimate_abundance.pl --transcripts /projectnb/davieslab/jfifer/Sym_Flow/ANNOTATIONS/Sym.fasta --est_method RSEM --aln_method bowtie2  --prep_reference --trinity_mode
#then run
/projectnb/davieslab/jfifer/Sym_Flow/trinityrnaseq-v2.9.0/util/align_and_estimate_abundance.pl --transcripts /projectnb/davieslab/jfifer/Sym_Flow/ANNOTATIONS/Sym.fasta --est_method RSEM --aln_method bowtie2  --trinity_mode --output_dir ./DE_output --samples_file /projectnb/davieslab/jfifer/Sym_Flow/COUNTS/samples.txt --seqType fq
