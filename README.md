# StaghornFlow
Scripts from GOING WITH THE FLOW: CORALS IN HIGH-FLOW ENVIRONMENTS CAN BEAT THE HEAT

# Transcriptome assembly, annotation and differential gene expression analysis pipeline
 These scripts were written by Dr. Bastian Bentlage and James Fifer and is a modified version of the pipeline available at https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Differential-Expression

# Transcriptome assembly
Step 1) Trim files using Trimmmomatic.sh script. Sequences were trimmed using TRIMMOMATIC (ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:35:-phred33)  (Bolger, Lohse, & Usadel, 2014), which removes low quality nucleotides with bp ≤35 and sequencing adapters, per default settings.

Step 2) Use tophat.sh on trimmed combined fq files. All quality-filtered paired reads were aligned against the publicly available A. digitifera genome (Shinzato et al., 2011), using the splice-junction mapper TopHat2 (Kim et al., 2013). Using the resulting BAM files, a reference transcriptome for A. cf. pulchra was assembled via the genome-guided version of the Trinity transcriptome assembler (Haas et al., 2013), with the A. digitifera genome (Shinzato et al., 2011) as a guide. 

Step 3) Use the parallelBlast.sh to split the fasta into one file per sequence in order to run in parallel. This script blasts against NCBI nt database and retains hits that match scleractinian. This creates a coralcontigs.txt file

Step 4) Filter out contigs from the coralcontigs.txt file with <300bp using the 300bpfilter.sh & contigs_to_fasta.py

Step 5) Filter out rRNA with rRNAblast.sh and remove.rRNA.sh.

# Annotation
Step 6) With this fasta file annotate using the script annotate.sh, blastx against the uniprot database with desired GO terms. Use script ParseAnnotations.py to only pull cnidarian GO terms from uniprot. Use FilterBLASTxResults.py to filter hits with e values < 1e-5 and Contigs_to_fasta.py to turn these results into fasta file. This creates the final reference transcriptome. 

# Differential gene expression analysis

Step 7) Format sample file according to samples.txt. Use the script DE.sh to create DESeq and edgeR directories. Use results.sh to find intersection between the two differential gene expression analyses and look for candidate genes. Additionally, to identify significant enrichment of GO terms for the up- and down-regulated genes we conducted a GO analysis using a Fisher’s exact test and rank-based adaptive clustering of GO terms (Wright, Aglyamova, Meyer, & Matz, 2015; https://github.com/z0on/GO_MWU). 


