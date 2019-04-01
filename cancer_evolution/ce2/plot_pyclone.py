#script for plotting data from a pyclone calculation
import matplotlib.pyplot as plt
import numpy as np
import sys
from plot_data import plot_data

dirname = sys.argv[1]
real = [10/12*100, 7/12*100, 4/12*100, 1/12*100] #reference data from assignment 2
cfs = []
sizes = []
f = open(dirname+'/tables/cluster.tsv', 'r') #file with results
f.readline()
for line in f:
	split = line.split()
	cfs.append(float(split[3])*100) #extract CCF
	sizes.append(float(split[2])) #extract #mutations
f.close()

sim = False
if len(sys.argv) > 2:
	if sys.argv[2] == 'sim': sim = True #simulated data
plot_data(cfs, sizes, dirname=dirname, sim=sim) #plot result

