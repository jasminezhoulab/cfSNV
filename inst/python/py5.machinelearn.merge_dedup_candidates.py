import sys
from operator import itemgetter

import pandas as pd
index_file = pd.read_csv(sys.argv[4], header=None, sep="\t")
chr_index = index_file[0]
chrom_dict = {}
for i in range(len(chr_index)):
    chrom_dict[chr_index[i]] = i


f = open(sys.argv[1], 'r')
occurred_pass = []
save = []
for line in f:
	sp = line.strip().split('\t')
	id = '-'.join(sp[0:4])
	if id not in occurred_pass:
		occurred_pass.append('-'.join(sp[0:4]))
		sp[1] = int(sp[1])
		sp.insert(1, sp[1]-1)
		sp.append("pass")
		sp.append(chrom_dict[sp[0]])
		save.append(tuple(sp))

f.close()

f = open(sys.argv[2], 'r')
occurred_hold = []
for line in f:
	sp = line.strip().split('\t')
	id = '-'.join(sp[0:4])
	if id not in occurred_pass and id not in occurred_hold:
		occurred_hold.append('-'.join(sp[0:4]))
		sp[1] = int(sp[1])
		sp.insert(1, sp[1]-1)
		sp.append("hold")
		sp.append(chrom_dict[sp[0]])
		save.append(tuple(sp))

f.close()


save.sort(key = itemgetter(2))
save.sort(key = itemgetter(18))

g = open(sys.argv[3], 'w')
for sp in save:
	sp = list(sp)
	sp[1] = str(sp[1])
	sp[2] = str(sp[2])
	g.write('\t'.join(sp[0:(len(sp)-1)]) + '\n')

g.close()

