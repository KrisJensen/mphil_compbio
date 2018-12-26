import numpy as np
import matplotlib.pyplot as plt
import sys

specs = ['H.sapiens', 'C.jacchus', 'D.rerio', 'R.norvegicus', 'M.musculus', 'G.gallus', 'F.catus', 'M.mulatta', 'I.punctatus', 'M.gallopavo']
N = len(specs)
name = 'RB1'
lsize = 14

if len(sys.argv) > 1:
	if sys.argv[1] == 'para':
		ori_specs = ['H.sapiens', 'M.musculus', 'M.mulatta', 'D.rerio']#', 'F.catus', 'M.mulatta', 'I.punctatus', 'M.gallopavo']
		isos = ['RB1', 'RBL1', 'RBL2']
		specs = []
		for spec in ori_specs:
			for iso in isos:
				specs.append(spec+'_'+iso)
		N = len(specs)
		name = 'para'



def get_sim(f):
	
	scores = []
	p0 = 'NA'

	for line in f:
		if line[:13] == 'global/global':
			return float(line.split()[4][:-1]), int(line.split()[3][:-1]) #percent sim, bit score
def get_dist(group1, group2, dists):
	distlist = []
	for s1 in group1:
		for s2 in group2:
			distlist.append( 0.50*(dists[s1][s2]+dists[s2][s1]))
	return np.mean(distlist)	

def cluster(specs, dists, name = 'RB1'):

	clusters = []
	groups = [ [spec] for spec in specs ]
	clusters.append(groups)
	connected = []
	c_dists = []

	while len(groups) > 1:
		smin = 0
		for i, g1 in enumerate(groups):
			for j, g2 in enumerate(groups):
				if j > i:
					s12 = get_dist(g1, g2, dists)
					if s12 > smin:
						smin = s12
						inds_min = (i, j)
		i, j = inds_min
		#print(groups[i], groups[j], smin)
		connected.append( [groups[i], groups[j]] )
		groups[i] = groups[i]+groups[j]
		groups.pop(j)
		clusters.append(groups)
		c_dists.append(100-smin)
		print(groups)

	xcoords = {}
	ycoords = {}
	for i, spec in enumerate(groups[0]):
		xcoords[spec] = i+1
		ycoords[spec] = 0
	scale = 0.8
	fig = plt.figure( figsize = (8*scale,6*scale) )
	ax = plt.axes()	
	plt.xlim(0.3, len(specs)+1)
	plt.ylim(0, 1.05*max(c_dists) )

	lwd = 4
	for i, d in enumerate(c_dists):
		x1, x2 = [ xcoords[connected[i][j][0]] for j in range(2) ]
		plt.plot([x1, x2], [d, d], 'k-', linewidth = lwd)
		plt.plot([x1, x1], [ycoords[connected[i][0][0]], d], '-k', linewidth = lwd)
		plt.plot([x2, x2], [ycoords[connected[i][1][0]], d], '-k', linewidth = lwd)

		ycoords[connected[i][0][0]] = d #keep track of how far up the tree we are at each point
		xcoords[connected[i][0][0]] = np.mean([x1, x2])
	for axis in ['right', 'top', 'bottom', 'left']:
		ax.spines[axis].set_visible(False)
	plt.ylabel('Divergence / %', fontweight='bold', fontsize = lsize)
	plt.yticks(fontweight='bold')
	ax.arrow(0.4, 0, 0, max(c_dists), width = 0.08, color = 'k', head_length = 0.4)
	plt.xticks(range(1, len(specs)+1), [ ' '.join(s.split('_')) for s in groups[0]],
		fontweight='bold', fontstyle = 'italic', rotation=45, ha = 'right')
	plt.tick_params(axis='both', labelsize = lsize)
	plt.savefig('/local/data/public/ktj21/GIA3/figures/'+name+'_dendrogram.png', dpi = 360, bbox_inches = 'tight')
	plt.show()
	
	return groups[0]

sim_mat = np.empty([N,N])
sims = {}

for i, spec1 in enumerate(specs):
	sims[spec1] = {}
	for j, spec2 in enumerate(specs):
		if spec1 != spec2:
			with open('blast/'+spec1+'_'+spec2+'.out', 'r') as f:
				sim, bit = get_sim(f)
				sim_mat[i,j] = sim
				print('\n'+' '.join(specs))
				print(np.round(sim_mat, 1))
				print(spec1, spec2, sim, bit)
				sims[spec1][spec2] = sim
#print(sims)
specs = cluster(specs, sims, name=name)
avg_mat = np.empty([N,N])
mi = np.inf
for i, spec1 in enumerate(specs):
	for j, spec2 in enumerate(specs):
		if spec1 == spec2: m = 100
		else:
			m = 0.5*(sims[spec1][spec2]+sims[spec2][spec1])
			if m < mi: mi = m
		avg_mat[i,j] = m
		avg_mat[j,i] = m
print('\n', np.round(avg_mat, 1))
#print(np.round(avg_norm, 3)) 	

plt.imshow(avg_mat, cmap='coolwarm', vmin = np.min(avg_mat), vmax = max(100, np.max(avg_mat)) )
plt.xticks(range(N), specs, rotation=40, ha='right', fontstyle = 'italic', fontweight='bold')
plt.yticks(range(N), specs, fontstyle = 'italic', fontweight='bold')
plt.tick_params(axis='both', labelsize = 15)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=15)
plt.savefig('/local/data/public/ktj21/GIA3/figures/'+name+'_percent_sim.png', bbox_inches = 'tight')
plt.show()

			
