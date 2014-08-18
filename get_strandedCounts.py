import sys, os

strand = dict([(x.strip().split()[0],x.strip().split()[1]) for x in open(sys.argv[1],'r')])

inFile = open(sys.argv[2])

orient = sys.argv[3]

for line in inFile:
        data = line.strip().split()

                if data[0] != 'no_feature':
                        tid = '.'.join(data[0].split('.')[:-1])
                        neg = data[1]
                        pos = inFile.next().strip().split()[1]

                        if strand[tid] == "+":
                                if orient == 'pos':
                                        print tid + '\t' + pos
                        else:
                                if orient == 'neg':
                                        print tid + '\t' + neg

                else:
                        break
