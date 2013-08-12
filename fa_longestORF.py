import sys
#requires BioPython
from Bio import SeqIO

#input is EMBOSS' getorf fasta output
inFile = open(sys.argv[1],'r')

gids = {}

for record in SeqIO.parse(inFile,'fasta'):
	gid = '_'.join(record.description.split()[0].split('_')[:-1])
	seq = str(record.seq)
	length = len(str(record.seq))

	start = record.description.split()[1][1:]
	end = record.description.split()[3][:-1]
	coord = start + "," + end

	if not gids.has_key(gid):
		gids[gid] = (seq, coord)
	else:
		if len(gids[gid][0]) < len(seq):
			gids[gid] = (seq, coord)

for gid, data in gids.items():
	print ">" + gid + "\t" + data[1]
	print data[0]
