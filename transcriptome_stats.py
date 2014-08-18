import sys,os
from Bio import SeqIO

inFile = open(sys.argv[1],'r')

lengths = []
for record in SeqIO.parse(inFile,'fasta'):
	lengths.append(len(str(record.seq)))

lengths.sort(reverse = True)

s = sum(lengths)
N50 = 0
count = 0
for l in lengths:
	count += l
	N50 = l
	if count > float(s) / 2:
		break

print 'number:', len(lengths)
print 'sum:', s
print 'mean:', (s / float(len(lengths)))
print 'N50:', l
print 'Top:'
print '\n'.join([str(x) for x in lengths[:10]])
