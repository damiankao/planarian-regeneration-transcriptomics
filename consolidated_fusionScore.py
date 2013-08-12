import sys, numpy, math, multiprocessing

#check if a position is within coordinate ranges
def isWithin(a, coord):
	if a >= coord[0] and a <= coord[1]:
		return True
	return False

#check if two cordinate ranges overlap
def overlap(coordA,coordB):
	if isWithin(coordA[0],coordB) or isWithin(coordA[1],coordB) or isWithin(coordB[0],coordA) or isWithin(coordB[1],coordA):
		return True
	return False

#convert coordinates into list of positions retaining "start" or "end" information
def sortCoord(d):
	coords = []
	for coord in d:
		coords.append(('s',coord[0]))
		coords.append(('e',coord[1]))

	coords.sort(key = lambda x : x[0], reverse = True)

	coords.sort(key = lambda x : x[1])

	return coords

#lexical parsing to find union of multiple coordinate ranges using sorted coordinates
def consensus(c):
	count = 0
	posA = 0
	out = []
	for pos in c:
		if count == 0:
			posA = pos[1]
		if pos[0] == 's':
			count += 1
		if pos[0] == 'e':
			count -=1

		if count == 0:
			out.append((posA, pos[1]))

	return out

#get total length given a list of coordinate range tuples
def getLength(coords):
	total = 0
	for coord in coords:
		total += coord[1] - coord[0] + 1
	return total

#given coordinates of alignments, give the longest alignment length
def getTopCoverage(consensusCoords, hsps):
	consensusLength = getLength(consensusCoords)

	hqueries = {}
	for hsp in hsps:
		if not hqueries.has_key(hsp[0]):
			hqueries[hsp[0]] = []
		hqueries[hsp[0]].append(hsp[1])

	queryLengths = []
	for query, coords in hqueries.items():
		consensusQuery = consensus(sortCoord(coords))
		queryLength = getLength(consensusQuery)
		queryLengths.append((query,queryLength,consensusQuery))

	queryLengths.sort(key = lambda x : x[1], reverse = True)

	return queryLengths[0]


queries = {}

#input file is a post-filtered tab delimited blast output with this format flag: -outfmt ‘6 qseqid qlen sseqid slen qstart qend sstart send evalue bitscore qframe’
inFile = open(sys.argv[1],'r')

currentQueries = {}
for line in inFile:
	data = line.strip().split()

	query = data[0]
	hit = data[2]
	qStart = int(data[4])
	qEnd = int(data[5])
	qLen = int(data[1])
	strand = '+'

	if data[-2][0] == "-":
		strand = "-"

	if not currentQueries.has_key(query):
		currentQueries[query] = []

	coord = [qStart,qEnd]
	coord.sort()

	currentQueries[query].append([hit, coord, float(data[9]), qLen, strand])

def multiThread(data, nprocs):
	def worker(subset, out_q):
		""" The worker function, invoked in a process. 'nums' is a
			list of numbers to factor. The results are placed in
			a dictionary that's pushed to a queue.
		"""
		outdict = {}
		for datum in subset:
			query = datum[0]
			hsps = datum[1]
			strands = []
			for item in hsps:
				strands.append(item[-1])

			strands = list(set(strands))
			consensusRegions = consensus(sortCoord([x[1] for x in hsps]))
			consensusLength = getLength(consensusRegions)
			topCoverage = getTopCoverage(consensusRegions, hsps)
			
			percentCoverage = float(topCoverage[1]) / consensusLength

			outdict[query]= (percentCoverage,consensusLength,hsps[0][-2],consensusRegions,topCoverage[2],strands)

		out_q.put(outdict)

	out_q = multiprocessing.Queue()
	chunksize = int(math.ceil(len(data) / float(nprocs)))
	procs = []

	for i in range(nprocs):
		p = multiprocessing.Process(
				target=worker,
				args=(data[chunksize * i:chunksize * (i + 1)],
					  out_q))
		procs.append(p)
		p.start()

	resultdict = {}
	for i in range(nprocs):
		resultdict.update(out_q.get())

	for p in procs:
		p.join()

	return resultdict

#multi-thread with 18 processes
queries = multiThread(currentQueries.items(),18)

print 'id\tpercentCov\ttopCovCoords\tconsensusCoord'
for query, data in queries.items():
	print query + "\t" + str(data[0]) + "\t" + ';'.join([str(x[0]) + "," + str(x[1]) for x in data[4]]) + "\t" + ';'.join([str(x[0]) + "," + str(x[1]) for x in data[3]]) + "\t" + ','.join(data[-1])
