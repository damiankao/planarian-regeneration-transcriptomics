import sys

#input is a tab delimited file where rows are transcripts and columns are expression value of various conditions. For example
#                head_0hours   head_6hours
#transcript01         100          200
#transcript02         400          100
#transcript03         900          100
inFile = open(sys.argv[1],'r')

libs = inFile.next().strip().split()[1:]

counts = {}

for lib in libs:
	counts[lib] = {}

for line in inFile:
	data = line.strip().split()
	tid = data[0]

	exps = data[1:]

	for i in range(len(libs)):
		lib = libs[i]
		exp = float(exps[i])

		counts[lib][tid] = exp

filtered = []
for lib in libs:
	vals = zip(counts[lib].keys(), counts[lib].values())
	total = sum([x[1] for x in vals])
	vals.sort(key = lambda x : x[1], reverse = True)

	for i in range(10):
		p = vals[i][1] / float(total) * 100
		if p >= 1:
			filtered.append(vals[i][0])

f = {}

for filterItem in filtered:
	if not f.has_key(filterItem):
		f[filterItem] = 1
	else:
		f[filterItem] += 1

for fItem, count in f.items():
	if count > 2:
		print fItem

