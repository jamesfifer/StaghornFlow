#first filter the blastoutput file so it only contains bp >300
awk '$4 >300 '{print $1} > blastoutput.300bp.plus.txt
#then use grep -F -f to filter out the IDs that are shared between
grep -F -f blastoutput.300bp.plus.txt CoralContigs.txt > filteredCoralContigs.txt
#then put the filteredtxt file through the Contigs_to_fasta.py script 
#then use rRNAblast.sh and remove.rRNA.sh to filter out RNA 

