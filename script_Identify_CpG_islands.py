#!/usr/bin/python

#script to identify CpG islands in each sequence of a given fasta using the method from https://www.bioinformatics.org/sms2/cpg_islands.html
#it works with 200 bp window size - default if nothing else given
#CpG Islands reports potential CpG island regions using the method described by Gardiner-Garden and Frommer (1987). The calculation is performed using a 200 bp window moving across the sequence at 1 bp intervals. CpG islands are defined as sequence ranges where the Obs/Exp value is greater than 0.6 and the GC content is greater than 50%. The expected number of CpG dimers in a window is calculated as the number of 'C's in the window multiplied by the number of 'G's in the window, divided by the window length. CpG islands are often found in the 5' regions of vertebrate genes, therefore this program can be used to highlight potential genes in genomic sequences.

#usage: python script_2023-12-11_identify_CpG-islands.py reference/GCA_900700415.2_Ch_v2.0.2_genomic.fa [optional: windowsize[default:200] stepsize[default:1]];

import sys
from Bio import SeqIO

#read in fa
fasta=SeqIO.parse(open(sys.argv[1]), 'fasta')

#read in optional parameters for window step size
try:
	win_step=int(sys.argv[3])
except IndexError:
	win_step=1

#ok now iterate through each fasta and determine CpG islands in a rolling window approach...
with open("CpGislands/"+sys.argv[1].split("/")[-1].split(".f")[0]+".CpGislands", 'w') as out:
	out.write("chr\tstart\tend\tC\tG\tGC\t%GC\tObs/Exp\n")
	for seq in fasta:

		#read in optional parameters for window start and end
		start=0
		try:
			win_size=int(sys.argv[2])
		except IndexError:
			win_size=200		

		while win_size<len(seq):
			window=seq.seq[start:win_size].upper()
			#print(window)
			Gcount=window.count("G")
			Ccount=window.count("C")
			GCcontent=round((window.count("G")+window.count("C"))/len(window)*100,2)
			GCcount=window.count("GC")
			try:
				obs_exp=round(GCcount/(Gcount*Ccount/len(window)),2)
			except ZeroDivisionError:
				obs_exp=0

			#write stats for current window to output if CpGisland (GC>50% and Obs/Exp > 0.6)
			if GCcontent>50 and obs_exp>0.6:
				out.write(seq.id+"\t"+str(start)+"\t"+str(win_size)+"\t"+str(Ccount)+"\t"+str(Gcount)+"\t"+str(GCcount)+"\t"+str(GCcontent)+"\t"+str(obs_exp)+"\n")			

			#prepare parameters for next window
			start=start+win_step
			win_size=win_size+win_step
