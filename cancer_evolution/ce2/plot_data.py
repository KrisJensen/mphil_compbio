#function for plotting clusters as barplotts

import matplotlib.pyplot as plt
import numpy as np
import sys

def plot_data(cfs, sizes, sim = False, real=[10/10*100, 7/10*100, 4/10*100, 1/10*100],
		dirname='.'):
	args = np.argsort(-np.array(cfs))
	cfs = np.array(cfs)[args] #CCFs
	sizes = np.array(sizes)[args] #number of mutations per cluster

	#print some summary data
	print('cfs:', cfs)
	print('clusters:', len(cfs))
	print('sizes:', sizes)
	S = np.sum(sizes)
	print('total sites:', S)
	print('>1:', np.sum(sizes > 1.5))
	dirname = dirname+'/'

	if sim: fsize = 7.5
	else: fsize = 12

	############plot inferred data
	ticks, va = [], []
	for i, f in enumerate(cfs): #need ticks for cellular fraction and CCF
		if sizes[i] > 1.5 or S < 30: ticks.append(str(np.round(f, 1)))
		else: ticks.append('')
		va.append(-0.03)
	fig = plt.figure(figsize = np.array([6,3])*0.8 )
	ax = plt.gca()
	plt.bar(cfs, sizes, width = 1.5)

	ax.set_xticks( cfs )
	ax.set_xticklabels( ticks, fontsize = fsize )
	plt.xlim(103,0)
	#ylims depend on the number of mutations considered
	if S > 150: smax = 126
	elif S > 50: smax = 70
	else: smax = 14
	if sim: smax = 20 
	plt.ylim(0, smax)

	plt.yticks(FontSize = 12)
	#plt.title('Column summary', FontSize=fs)
	plt.xlabel('Cancer Cell Fraction (%)', FontSize=14)
	plt.ylabel('# mutations', FontSize=14)
	plt.savefig(dirname+'cols.png', bbox_inches = 'tight', dpi=360)
	plt.close()

	#############plot reference data

	if sim: #if working with simulated data, reference data is from assingment 1
		real = [100, 75, 50, 30, 20, 15] 
		realsizes = [5, 11, 3, 6, 6, 9]
	else: realsizes = []
	ticks, va = [], []
	for i, f in enumerate(real): #need ticks for cellular fraction and CCF
		ticks.append(str(np.round(f, 1)))
		va.append(-0.03)
		if not sim: realsizes.append(180/4)
	fig = plt.figure(figsize = np.array([6,3])*0.8 )
	ax = plt.gca()
	plt.bar(real, realsizes, width = 1.5)

	ax.set_xticks( real )
	ax.set_xticklabels( ticks, fontsize = fsize )
	plt.xlim(103,0)
	plt.ylim(0, smax)
	plt.yticks(FontSize = 12)
	plt.xlabel('Cancer Cell Fraction (%)', FontSize=14)
	plt.ylabel('# mutations', FontSize=14)
	plt.savefig(dirname+'realcols.png', bbox_inches = 'tight', dpi=360)
	plt.close()


