import sys
W = int(sys.argv[1])
fname = sys.argv[2]
gname = sys.argv[3]
f = open(fname, 'r')
g = open(gname, 'w')
for line in f:
	sp = line.strip().split('\t')
	sp[1] = str(int(sp[1])-W)
	sp[2] = str(int(sp[2])+W)
	g.write('\t'.join(sp[0:3] + [sp[4]]) + '\n')

g.close()
f.close()
