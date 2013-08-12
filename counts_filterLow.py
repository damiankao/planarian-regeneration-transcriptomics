import sys, numpy

#input is a tab delimited file where rows are transcripts and columns are expression value of various conditions. For example
#                head_0hours   head_6hours
#transcript01         100          200
#transcript02         400          100
#transcript03         900          100

inFile = open(sys.argv[1],'r')

print inFile.next().strip()

for line in inFile:
	data = line.strip().split()

	vals = [float(x) for x in data[1:]]	
	check = False

	for val in vals:
		if val >= 20:
			check = True
	if check:
		print line.strip()