import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os, glob

dfiles = sorted(glob.glob('*.txt'))

SMALL_SIZE = 12
MEDIUM_SIZE = 14
BIGGER_SIZE = 16

#plt.style.use('bmh')
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
plt.rc('figure', figsize=(8,8))

def splot():
	try:
		plotnumber = int(input('Enter plot number: '))
		pt = pd.read_csv(dfiles[plotnumber], sep='\t', header=None, skiprows=1, engine='python')
		for n in range(0,len(pt.iloc[0,:])//2):
			plt.plot(pt.iloc[:,2*n], pt.iloc[:,2*n+1], label=n)
		plt.axhline(0, 0, 1, c = 'black')
		plt.xlabel('$E$/V vs Pd/H')
		plt.ylabel(r'$i\,/\,\mathrm{\mu A}$')
		plt.legend()
		plt.show()
	except ValueError:
		print('Enter an integer')
		splot()	

def exp_plots():
	options = {'n':'n', 'N':'n', 'y':'y', 'Y':'y'}
	selection = input("Save all plots? (y/N) ")
	ret = options.get(selection)
	if not ret:
		return exp_plots()
	if ret == 'y':
		for k in range(0,len(dfiles)):
			pt = pd.read_csv(dfiles[k], sep='\t', header=None, skiprows=1, engine='python')
			for n in range(0,len(pt.iloc[0,:])//2):
				plt.plot(pt.iloc[:,2*n], pt.iloc[:,2*n+1], label=n)
			plt.axhline(0, 0, 1, c = 'black')
			plt.xlabel('$E$/V vs Pd/H')
			plt.ylabel(r'$i\,/\,\mathrm{\mu A}$')
			plt.legend()
			plt.savefig('voltagramas%02d.svg' % (k))
			plt.clf()
	else:
		return print('...plots not saved')

choic = input("Choose a plot? (y): ")
if choic == 'y' :
	splot()

exp_plots()
