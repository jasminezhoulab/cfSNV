import sys
import pandas as pd
x = pd.read_csv(sys.argv[1], header = None, dtype = {0: str}, sep='\t')

x[11] = x[1]-4
x[12] = x[1]+5

from collections import defaultdict
INDEL_COUNT = 3
INDEL_VAF = 0.02
indel_bins = defaultdict(list)
for i in range(x.shape[0]):
    if x[5][i] > INDEL_COUNT or x[7][i]+x[9][i] > INDEL_COUNT or (float(x[6][i]) > 0 and float(x[5][i])/float(x[6][i]) > INDEL_VAF) or (float(x[8][i]) + float(x[10][i]) > 0 and float(x[7][i]+x[9][i])/float(x[8][i]+x[10][i]) > INDEL_VAF):
        current_start = x[11][i]
        current_end = x[12][i]
        current_type = x[2][i]
        current_chrom = x[0][i]
        indel_bins[current_chrom].append(((current_start, current_end), current_type))
    

indel_merged_bins = []

index_file = pd.read_csv(sys.argv[3], header=None, sep="\t")
chrom_list = index_file[0]

for i in chrom_list:
    flag = 0
    if i not in indel_bins.keys():
        continue
    for j in indel_bins[i]:
        if flag == 0:
            current_type = j[1]
            current_region = j[0]
            flag = 1
        else:
            if j[1] != current_type:
                indel_merged_bins.append((i, current_region, current_type))
                current_type = j[1]
                current_region = j[0]
            elif current_region[0] < j[0][0] and current_region[1] >= j[0][0]:
            #   print current_region, j[0]
                current_region = (min(current_region[0], j[0][0]), max(current_region[1], j[0][1]))
            else:
                indel_merged_bins.append((i, current_region, current_type))
                current_type = j[1]
                current_region = j[0]
    indel_merged_bins.append((i, current_region, current_type))



output = open(sys.argv[2], 'w')

for i in indel_merged_bins:
    output.write(i[0] + '\t' + str(i[1][0]) + '\t' + str(i[1][1]) + '\t' + i[2] + '\n')

output.close()
