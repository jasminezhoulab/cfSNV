import sys

overlap = open(sys.argv[2], 'w')
nonoverlap = open(sys.argv[3], 'w')
doverlap = {}
dnonoverlap = {}
noverlap = 0
nnonoverlap = 0
varoverlap = ""
startoverlap = 0
varnonoverlap = ""
startnonoverlap = 0
with open(sys.argv[1], 'r') as raw:
	for line in raw:
		sp = line.strip().split('\t')
		var = '\t'.join(sp[61:64])
		if sp[36] == "NA":
			if varnonoverlap == var:
				nonoverlap.write(line)
				nnonoverlap += 1
			else:
				dnonoverlap[varnonoverlap] = [str(startnonoverlap), str(nnonoverlap)]
				varnonoverlap = var
				nonoverlap.write(line)
				startnonoverlap = nnonoverlap
				nnonoverlap += 1
		if sp[36] != "NA":
			if varoverlap == var:
				overlap.write(line)
				noverlap += 1
			else:
				doverlap[varoverlap] = [str(startoverlap), str(noverlap)]
				varoverlap = var
				overlap.write(line)
				startoverlap = noverlap
				noverlap += 1

overlap.close()
nonoverlap.close()
#print(noverlap, nnonoverlap)

#overlap = open(sys.argv[4], 'w')
#nonoverlap = open(sys.argv[5], 'w')
#for i in dnonoverlap.keys():
#	if i == "":
#		continue
#	nonoverlap.write(i + '\t' + '\t'.join(dnonoverlap[i]) + '\n')

#nonoverlap.close()
#for i in doverlap.keys():
#	if i == "":
#		continue
#	overlap.write(i + '\t' + '\t'.join(doverlap[i]) + '\n')

#overlap.close()

