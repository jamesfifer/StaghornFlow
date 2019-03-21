


#for sample in `ls ./*R1_Combined.fq`
#do
#	outdir=`echo ${sample/_R1_Combined.fq/}`
#	nice -n 6 tophat -o $outdir /media/RAID/james/Staghorn_GeneExpression/In_Situ/DigitiferaRefGenome/Adig $sample ${sample/R1/R2} &
#done
#Takes all files in all downhill directories with name and copies with unique names to current directory
#find ./ -name 'accepted_hits.bam' -exec cp --backup=numbered -t ./ {} + 
#samtools merge merged.bam *bam*

#Trinity recommends to make sure the file is coordinate sorted by running 'samtools sort' on it

#Run genome guided trinity assembly
Trinity --genome_guided_bam /media/RAID/james/Staghorn_GeneExpression/Ref_Transcriptome_All/tophat_out/accepted_hits.bam \
       --CPU 2 --max_memory 1G --genome_guided_max_intron 5000
