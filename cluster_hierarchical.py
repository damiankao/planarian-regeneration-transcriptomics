import sys, numpy, scipy
import scipy.cluster.hierarchy as hier
import scipy.spatial.distance as dist
import pylab

order = '''0 hr head
6 hr head
12 hr head
24 hr head
36 hr head
48 hr head
72 hr head
0 hr tail
6 hr tail
12 hr tail
24 hr tail
36 hr tail
48 hr tail
72 hr tail'''.strip().split('\n')

#input is a tab delimited m x n matrix of observations where rows are genes and columns are conditions. Row and column headers are included in the input file.
inFile = open(sys.argv[1],'r')
colHeaders = inFile.next().strip().split()[1:]
rowHeaders = []
dataMatrix = []
for line in inFile:
	data = line.strip().split()
	rowHeaders.append(data[0])
	dataMatrix.append([float(x) for x in data[1:]])

dataMatrix = numpy.array(dataMatrix)

print 'calculate correlation matrix'
distanceMatrix = dist.pdist(dataMatrix,'correlation')

print 'calculate linkage matrix'
linkageMatrix = hier.linkage(distanceMatrix,method='complete')

def lablef(i):
	return order[i]

fig = pylab.figure(figsize=(10,7))
Z1 = hier.dendrogram(linkageMatrix,leaf_label_func=lablef,orientation='left')
fig.savefig('dendrogram.pdf')