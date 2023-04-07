import sys
import numpy as np

def depth_range(region_depth, output):
  with open(region_depth, "r") as f:
    depths = f.readlines()
    depths = [float(i.strip()) for i in depths]
    avr = np.average(depths)
    qtl = np.quantile(depths, [0.5, 0.95])
    
  output_f = open(output, 'w')
  output_f.write(str(avr) + '\n')
  output_f.write(str(qtl[0]) + '\n')
  output_f.write(str(qtl[1]) + '\n')
  output_f.close()
    
    # return avr, qtl[0], qtl[1]

region_depth = sys.argv[1]
region_depth_output = sys.argv[2]
depth_range(region_depth, region_depth_output)
