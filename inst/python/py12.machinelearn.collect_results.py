import sys

path = sys.argv[1]
plasma = sys.argv[2]
sample = path + plasma

#threshold = [0.8, 0.85, 0.9, 0.95] # thresholds for RF20

threshold = [0.2, 0.3, 0.4, 0.5] # thresholds for RF41


f = open(sample+".bed")
variant_list = {}
for line in f:
	sp = line.strip().split('\t')
	variant_list['-'.join([sp[0], sp[2]])] = [sp[4], int(sp[13]), 0, 0] + [0]*(2*len(threshold))

f.close()


g = open(sample+".paired-reads.qsort.nonoverlap.features")
h = open(sample+".paired-reads.qsort.nonoverlap.features.expand.RFpred.csv")
for line in g:
	sp = line.strip().split('\t')
	var = sp[len(sp)-6] + '-' + sp[len(sp)-5]
	pred = float(h.readline().strip())
	variant_list[var][2] += 1
	for i in range(len(threshold)):
		if pred > threshold[i]:
			variant_list[var][4+i*2] += 1

g.close()
h.close()

g = open(sample+".paired-reads.qsort.overlap.features")
h = open(sample+".paired-reads.qsort.overlap.features.expand.RFpred.csv")
for line in g:
	sp = line.strip().split('\t')
	var = sp[len(sp)-6] + '-' + sp[len(sp)-5]
	pred = float(h.readline().strip())
	variant_list[var][3] += 1
	for i in range(len(threshold)):
		if pred > threshold[i]:
			variant_list[var][5+i*2] += 1

g.close()
h.close()


x = open(sample+".after_machine_learn", 'w')
f = open(sample+".bed")
for line in f:
	sp = line.strip().split('\t')
	#sp[1] = str(int(sp[1]) + 1)
	var = '-'.join([sp[0], sp[2]])
	out = '\t'.join(sp) + '\t' + '\t'.join([str(i) for i in variant_list[var][2:]]) + '\n'
	x.write(out)

f.close()
x.close()


