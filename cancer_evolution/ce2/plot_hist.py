import matplotlib.pyplot as plt
import numpy as np

vals = []
f = open('PyClone_sim.tsv', 'r')
f.readline()
for line in f:
	split = line.split()
	vals.append(float(split[2])/(float(split[2])+float(split[1])))

f.close()

plt.hist(vals)
plt.xlabel('$f_B$')
plt.ylabel('frequency')
plt.savefig('vafhist.png')
plt.show()


