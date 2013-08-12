import sys

#tab delimited file of nucleotide lengths where column 1 is transcript id and column 2 is length
nucLengths = dict([(x.strip().split()[0],int(x.strip().split()[1])) for x in open(sys.argv[1],'r')])
#tab delimited file of longest ORF lengths where column 1 is transcript id and column 2 is length
orfLengths = dict([(x.strip().split()[0],int(x.strip().split()[1])) for x in open(sys.argv[2],'r')])
#clustering file where each line is a cluster of comma separated member transcript ids
clusterFile = open(sys.argv[3],'r')

for line in clusterFile:
	members = line.strip().split(',')
	if len(members) == 1:
		print members[0]
	else:
		data = []
		for member in members:
			if orfLengths.has_key(member):
				data.append((nucLengths[member],orfLengths[member], member))

		data.sort(key = lambda x : x[0], reverse = True)

		nuc_longest = data[0][0]

		nucs = []
		for item in data:
			#pool sequences that are at least 90% of the CAP3 contig length
			if item[0] / float(nuc_longest) >= 0.9:
				nucs.append(item)

		nucs.sort(key = lambda x : x[1], reverse = True)

		print nucs[0][-1]

