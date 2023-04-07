import sys
import numpy as np

def estimate_sequencing_depth(target_bed, base_mapped, output):
	target_bed = open(target_bed)
	base_mapped = open(base_mapped)
	n_read_base = float(base_mapped.readline().strip())
	n_region_base = 0
	for line in target_bed:
		sp = line.strip().split('\t')
		n_region_base = n_region_base + int(sp[2]) - int(sp[1])

	base_mapped.close()
	target_bed.close()
	
	depth = round(float(n_read_base)/float(n_region_base) * 0.8, 9)
		
	output = open(output, 'w')
	output.write(str(depth) + '\n')
	output.close()
	
	# return depth
    
target_bed = sys.argv[1]
base_mapped = sys.argv[2]
output = sys.argv[3]

estimate_sequencing_depth(target_bed, base_mapped, output)
