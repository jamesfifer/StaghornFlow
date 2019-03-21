#!/bin/bash
#bowtie2-build Sym_COIDB.final.fasta Sym_COIDB

#/usr/local/src/cdhit/cd-hit-est -i Sym_COIDB2.fasta -c 1 -o Sym_COIDB.final.fasta


for sample in `ls /media/RAID/james/Staghorn_GeneExpression/Ex_Situ/*R1_Combined.fq`
do
dir="/media/RAID/james/Staghorn_GeneExpression/Ex_Situ/"
base=$(basename $sample "R1_Combined.fq")
bowtie2 -x ./Sym_COIDB -1 ${dir}/${base}R1_Combined.fq -2 ${dir}/${base}R2_Combined.fq --no-unal -S ~/Staghorn_GeneExpression/Ex_Situ/Sym_Cladistics/${base}.COI.sam
done

for sample in `ls /media/RAID/james/Staghorn_GeneExpression/Ex_Situ/Sym_Cladistics/*.sam`
do
dir="/media/RAID/james/Staghorn_GeneExpression/Ex_Situ/Sym_Cladistics/"
base=$(basename $sample ".sam")
samtools view -bS ${dir}/${base}.sam >${dir}/${base}.bam
done 



for sample in `ls /media/RAID/james/Staghorn_GeneExpression/Ex_Situ/Sym_Cladistics/*.bam`
do
dir="/media/RAID/james/Staghorn_GeneExpression/Ex_Situ/Sym_Cladistics/"
base=$(basename $sample ".bam")
samtools bam2fq ${dir}/${base}.bam > ${dir}/${base}.fq 
done

for sample in `ls /media/RAID/james/Staghorn_GeneExpression/Ex_Situ/Sym_Cladistics/*.COI.fq`
do
dir="/media/RAID/james/Staghorn_GeneExpression/Ex_Situ/Sym_Cladistics/"
base=$(basename $sample ".COI.fq")
nice -n 6 Trinity --seqType fq --single ${dir}/${base}.COI.fq --max_memory 100G --CPU 40 --output ${dir}/${base}.trinityout
done



