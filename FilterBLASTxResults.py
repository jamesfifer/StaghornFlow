#!/usr/bin/env python

with open('BlastxAnnotated.GO.txt', 'r') as BLAST:
	for Line in BLAST:
		if float(Line.split()[3]) <=  1e-5:
			print Line.rstrip()
