import sys
#requires BioPython
from Bio import SeqIO

#input file is a fasta file
inFile = open(sys.argv[1],'r')
for record in SeqIO.parse(inFile,'fasta'):
	print record.id + "\t" + str(len(str(record.seq)))
