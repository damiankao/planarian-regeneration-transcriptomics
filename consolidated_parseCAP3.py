import sys

#input is the cap3 stdout
inFile = open(sys.argv[1],'r')
inFile.next()
inFile.next()
inFile.next()
inFile.next()
inFile.next()

header = "Contig" + inFile.next().strip().split()[2]
block = []

clusters = {}

inCluster = {}
repCluster = {}
for line in inFile:
	if line.find('DETAILED DISPLAY OF CONTIGS') != -1:
		clusters[header] = []
		for bline in block[:-1]:
			clusters[header].append(bline.strip().split()[0])
		break
	else:
		if line[0] == "*":
			clusters[header] = []
			for bline in block:
				clusters[header].append(bline.strip().split()[0])
			header = "Contig" + line.strip().split()[2]
			block = []
		else:
			block.append(line.strip())

for rep, members in clusters.items():
	print rep + "\t" + ','.join(members)
