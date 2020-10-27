#THIS SPLITS A LARGE  FILE INTO 1 FILE PER SEQUENCE
# awk '/^>/ {if(x>0) close(outname); x++; outname=sprintf("_%d.fa",x); print > outname;next;} {if(x>0) print >> outname;}' *.fasta


#this performsmany blasts at once. USE:
#find -type f -name '*.fa' | parallel -a - blastn -query {} -db /media/RAID/DATABASES/NCBI/nt -evalue 1e-5 -outfmt 6 -max_target_seqs 1 -max_hsps 1 -out {.}.out

# Parse the NCBI Taxonomy to retain only those hits that match scleractinians
for i in ./*out
do
	if [ -s "$i" ]; then
		Contig_ID=`awk '{ print $1 }' $i`
		NCBI_Acc=`awk '{ print $2 }' $i`
		Taxon=`esearch -db nuccore -query "$NCBI_Acc" | elink -target taxonomy | efetch -format uid`
		Taxonomy=`/usr/local/src/anaconda_ete/python-scripts/ete3  ncbiquery --search "$Taxon" --info | grep -v '#' | awk '{ print $6 }'`
		if [[ "$Taxonomy" = *"Scleractinia"* ]]; then
			echo "$Contig_ID" >> "CoralContigs.txt"
		fi
	fi
done

#for symbiont 
#blast search with these parameters (did nr here because of lower number of hits to nt)
#find -type f -name '*.fa' | parallel -a - blastx -query {} -db /media/RAID/DATABASES/NCBI/nr -evalue 1e-5 -outfmt 6 -max_target_seqs 1 -max_hsps 1 -out {.}.out

#remove .out files that did not return hits
# find . -name '*.out' -size 0 -print0 | xargs -0 rm

# Parse the NCBI Taxonomy to retain only those hits that match Symbiodiniaceae
for i in ./*out
do
       if [ -s "$i" ]; then
               Contig_ID=`awk '{ print $1 }' $i`
               NCBI_Acc=`awk -F '|' '{ print $4 }' $i`
               Taxon=`esearch -db protein -query "$NCBI_Acc" | elink -target taxonomy | efetch -format uid`
#               Taxon=`esearch -db nuccore -query "$NCBI_Acc" | elink -target taxonomy | efetch -format uid`
               Taxonomy=`/usr/local/src/anaconda_ete/python-scripts/ete3  ncbiquery --search "$Taxon" --info | grep -v '#'`
               if [[ "$Taxonomy" = *"Dinophyceae"* ]]; then
                       echo "$Contig_ID" >> "Symcontigs.txt"
                       echo -e "$Contig_ID" "\n"
               fi
       fi
done
