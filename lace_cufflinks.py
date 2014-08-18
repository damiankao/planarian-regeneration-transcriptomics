import sys
from collections import defaultdict

inFile = open(sys.argv[1],'r')

refs = defaultdict(list)
for line in inFile:
	data = line.strip().split('\t')
	if data[2] == 'exon':
		refs[data[0]].append((int(data[3]), int(data[4])))

count = 0
for ref, coords in refs.items():
	coords.sort(key = lambda x : x[0])
	for coord in coords:
		print '\t'.join([ref,'bioscope','exon',str(coord[0]),str(coord[1]),'.','+','.','gene_id "bioscope.pos.' + str(count) + '"; transcript_id "bioscope.pos' + str(count)'";'])

	for coord in coords:
		print '\t'.join([ref,'bioscope','exon',str(coord[0]),str(coord[1]),'.','-','.','gene_id "bioscope.neg.' + str(count) + '"; transcript_id "bioscope.neg' + str(count)'";'])

	count += 1
