from parameter import *
from _jenks import *
from _filter import *
import sys
from collections import defaultdict
from operator import itemgetter
from scipy.stats import binom

output_dir = sys.argv[1]
plasma = sys.argv[2]
results_dir = sys.argv[3]
MIN_HOLD_SUPPORT_COUNT = int(sys.argv[5])
MIN_PASS_SUPPORT_COUNT = int(sys.argv[6])

sample = output_dir + plasma

VAF_list = []
f = open(sample + ".after_cluster", 'r')
for line in f:
	sp = line.strip().split('\t')
	support_cnt = int(sp[26]) + int(sp[27])
	if sp[17] == "hold" and support_cnt >= MIN_HOLD_SUPPORT_COUNT:
		if int(sp[7]) > 50 and int(sp[12]) > 50:
			VAF_list.append(float(sp[15]))
	if sp[17] == "pass" and support_cnt >= MIN_PASS_SUPPORT_COUNT:
		if int(sp[7]) > 50 and int(sp[12]) > 50:
			VAF_list.append(float(sp[15]))

f.close()

jenks_estimate, include_number, include_group = final_estimation_with_Jenks(VAF_list)

dict = defaultdict(list)
VAF_list = []
f = open(sample + ".after_cluster", 'r')
a = 0
b = 0
c = 0
d = 0
for line in f:
	sp = line.strip().split('\t')
	support_cnt = int(sp[26]) + int(sp[27])
	total_cnt = int(sp[12])
	vaf = float(sp[15])
	if sp[17] == "hold" and  support_cnt > MIN_HOLD_SUPPORT_COUNT:
		output_seq = sp[0]+'\t' + sp[2] + '\t' + '.' + '\t' + sp[3] + '\t' + sp[4] + '\t' + sp[11] + '\t' + 'PASS' + '\t' + 'VAF=' + sp[10] + '\n'
		dict[sp[0]].append( ( int(sp[2]), output_seq ) )
		a += 1
		if int(sp[7]) > 50 and int(sp[12]) > 50:
			VAF_list.append(float(sp[10]))
	elif sp[17] == "hold" and support_cnt > 5 and  binom.cdf(support_cnt, total_cnt, jenks_estimate/2.0) > 0.4:
		output_seq = sp[0]+'\t' + sp[2] + '\t' + '.' + '\t' + sp[3] + '\t' + sp[4] + '\t' + sp[11] + '\t' + 'PASS' + '\t' + 'VAF=' + sp[10] + '\n'
		dict[sp[0]].append( ( int(sp[2]), output_seq ) )
		b += 1
		if int(sp[7]) > 50 and int(sp[12]) > 50:
			VAF_list.append(float(sp[10]))
	if sp[17] == "pass" and support_cnt > MIN_PASS_SUPPORT_COUNT:
		output_seq = sp[0]+'\t' + sp[2] + '\t' + '.' + '\t' + sp[3] + '\t' + sp[4] + '\t' + sp[11] + '\t' + 'PASS' + '\t' + 'VAF=' + sp[10] + '\n'
		dict[sp[0]].append( ( int(sp[2]), output_seq ) )
		c += 1
		if int(sp[7]) > 50 and int(sp[12]) > 50:
			VAF_list.append(float(sp[10]))
	elif sp[17] == "pass" and support_cnt > 3 and binom.cdf(support_cnt, total_cnt, jenks_estimate/2.0) > 0.4:
		output_seq = sp[0]+'\t' + sp[2] + '\t' + '.' + '\t' + sp[3] + '\t' + sp[4] + '\t' + sp[11] + '\t' + 'PASS' + '\t' + 'VAF=' + sp[10] + '\n'
		dict[sp[0]].append( ( int(sp[2]), output_seq ) )
		d += 1
		if int(sp[7]) > 50 and int(sp[12]) > 50:
			VAF_list.append(float(sp[10]))

f.close()

#print(a, b, c, d)



output = open(results_dir + "/" + plasma + ".variant_report.txt", 'w')
#output.write("CHROM\tPOSITION\tID\tREF\tVAR\tSCORE\tVAF\n")

import pandas as pd
index_file = pd.read_csv(sys.argv[4], header=None, sep="\t")
chr_index = index_file[0]
#for i in ["chr" + str(i) for i in range(1,23)] + ["chrX", "chrY"]:
for i in chr_index:
	current = sorted(dict[i], key=itemgetter(0))
	for item in current:
		output.write(item[1])

output.close()

jenks_estimate, include_number, include_group = final_estimation_with_Jenks(VAF_list)
with open(results_dir + "/" + plasma + ".jenks_estimate", 'w') as TF:
  TF.write(str(round(jenks_estimate, 12)) + '\t' + str(include_number) + '\t' + str(include_group) + '\n')
