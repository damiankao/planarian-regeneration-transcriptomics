import sys

#each line of the input file is a cluster of transcripts represented by a comma separated list of member ids.
inFile = open(sys.argv[1],'r')

for line in inFile:
	#each transcript id is in the format transcript source and an incremeted number. For example: aboobaker.1, bartscherer.123. The source can be obtained by id.split('.')[0]
	sources = set([x.split('.')[0] for x in line.strip().split(',')])
	if len(sources) > 1:
		print line.strip()
