#! /usr/bin/python

#Input: FASTA file name
# Searches a FASTA file for a pattern (from a file with an ID on each line) and prints the corresponding fasta entry to a new file. 


import sys, Bio


from Bio import SeqIO



fastaFile = sys.argv[1]
ids_file = sys.argv[2]
newfile = sys.argv[3]

# Put all Ensembl IDs in a set
IDset = set(line.strip() for line in open(ids_file))

HitFastaRecord = []

from Bio import SeqIO
# for each fasta record in the fasta file
for seq_record in SeqIO.parse(fastaFile, "fasta"):
	seq_recordID = seq_record.description.split()[0]
	if seq_recordID in IDset:
		HitFastaRecord.append(seq_record)
		print "wrote fasta record for:" + seq_recordID
	else:
		continue



			
output_file = open(newfile, "w")
SeqIO.write(HitFastaRecord, newfile, "fasta")
output_file.close()





