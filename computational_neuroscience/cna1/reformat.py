
from subprocess import call

fnames = ['firing_rates', 'firing_rates_finegrained']

for name in ['V', 'n', 'm', 'h', 'gating']:
	fnames+=['hh_'+name, 'hh_jump_'+name]

for name in fnames:
	call(['sips -s format png '+name+'.pdf --out '+name+'.png'], shell=True)


