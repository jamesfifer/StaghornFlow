
#module load boost/1.58.0
#module load bowtie2
#module load tophat/2.1.1
#module load samtools/1.9
#Below uses the digitifera genome as an example, but for the symbiont the same thing was done (subbed goreaui for Adig genome)

#use bowtie2 to build the databases (sym and host separately) before running this script


#for sample in `ls ./*R1_Combined.fq`
#do
#	outdir=`echo ${sample/_R1_Combined.fq/}`
#	nice -n 6 tophat -o $outdir /media/RAID/james/Staghorn_GeneExpression/In_Situ/DigitiferaRefGenome/Adig $sample ${sample/R1/R2} &
#done


#Takes all files in all downhill directories with name and copies with unique names to current directory
#find ./ -name 'accepted_hits.bam' -exec cp --backup=numbered -t ./ {} + 

#samtools merge merged.bam *bam*

#Trinity recommends to make sure the file is coordinate sorted by running 'samtools sort' on it
#samtools sort merged.bam > accepted_hits.merged.sorted.bam

#Run genome guided trinity assembly
module load java/1.8.0_181
module load bowtie2
module load samtools
module load jellyfish
module load salmon
module load python3/3.6.5
module load trinity/2.8.4


Trinity --genome_guided_bam /media/RAID/james/Staghorn_GeneExpression/Ref_Transcriptome_All/tophat_out/accepted_hits.bam \
       --CPU 2 --max_memory 1G --genome_guided_max_intron 5000
