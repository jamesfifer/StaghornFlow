#!/usr/bin/env python
#Make sure you have the uniprot database created from uniprot fasta, plus file with all go terms (.tab), and the blast outputs (blastx and blastp)

import os, sys

# Ceate a dictionary that contains all protein annotations, including
# GO terms from Uniprot
AnnotationDict = {}
with open('uniprot-cnidaria.GO.tab', 'r') as Annotations:
	for Annotation in Annotations:
		if 'Entry' not in Annotation:
			UniprotID = Annotation.rstrip().split('\t')[1]
			Annotation = Annotation.rstrip().split('\t')[2:]
			AnnotationDict[UniprotID] = Annotation

# Parse the BLAST results and print an annotation report that links
# BLAST results to the Uniprot annotations
with open('blastx.outfmt6', 'r') as Blast:
	for Line in Blast:
		BlastResults = Line.rstrip().split()
		UniprotID = BlastResults[1].split('|')[-1]
		Annotation = ('\t').join(AnnotationDict[UniprotID])
		sys.stdout.write('%s\t%s\t%s\t%s\t%s\n' % (BlastResults[0], UniprotID, BlastResults[2], BlastResults[-2], Annotation))
