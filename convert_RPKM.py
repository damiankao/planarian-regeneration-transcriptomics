import sys
from Bio import SeqIO

faFile = open(sys.argv[1],'r')

lengths = {}

for record in SeqIO.parse(faFile,'fasta'):
	lengths[record.id] = len(str(record.seq))

countFile = open(sys.argv[2],'r')
tids = {}
for line in countFile:
	data = line.strip().split('\t')
	tids[data[0]] = int(data[1])

libSize = float(sum([x[1] for x in tids]))

for tid, count in tids.items();
	rpkm = count / libSize / lengths[tid] * 1000
	print tid + '\t' + str(rpkm)
