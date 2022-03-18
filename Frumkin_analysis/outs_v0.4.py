'''
Function that fits a number of results of electrochemical impedance
spectroscopy and prints the fitted parameters as output.
'''
import numpy as np
import os, glob, levmpy, sys
import pandas as pd
import matplotlib.pyplot as plt

#%02d stands for a number of two digits and the % sign indicates the source
#getcwd prints the current directory path
#Suppress verbose fit info

fp = pd.DataFrame() #define empty dataframe for parameters
fd = pd.DataFrame() #frequency
zexp = pd.DataFrame() #initial impedance data
yeee = pd.DataFrame() #initial admittance data

def fit_eqs():
	options = {'z':'z', 'Z':'Z', 'y':'y', 'Y':'Y'} #dictionary to select mode
	selection = input("Choose Z for impedance input or Y for admittance input ")
	ret = options.get(selection) #get dictionary
	if not ret: #repeat function if incorrect commands are given
		return fit_eqs()
	if ret == 'z' or ret == 'Z':
		print('Impedance')
		ifiles = sorted(glob.glob('*1.txt')) #reads all files ending in .txt
		for k in range(0,len(ifiles)-1): #select INPUT files
			expt = levmpy.Experiment(ifiles[k+1], outsuffix="OUTZ%02d" % (k+1), path=os.getcwd()) #need to add an if statement to choose if outfiles are wanted
			expt.fit() #fit command
			fp.loc[:,k] = expt.x #creates fit parameters database
			fd.loc[:,k] = expt.freq #creates frequency database
			zexp.loc[:,k] = expt.y1 + expt.y2*1j #experimental impedance data
		fp.to_csv('FitParametersZZRRF.csv', index=False) #exports fit parameters database to CSV file without indexes
	if ret == 'y' or ret == 'Y':
		print('Admittance')
		ifiles = sorted(glob.glob('*1.txt'))
		for k in range(0,len(ifiles)-1): #cycle that changes input files names and makes the fit in admittance
			with open(ifiles[k+1], "rt") as fin:
				with open(ifiles[k+1][:ifiles[k+1].find('.')] + 'z.tmp', "wt") as fout: #select strings before the dot to establish name
					for line in fin:
						fout.write(line.replace('ZZRRF','ZYRRF')) #Z: input, Y: fit, RR: Nyquist result, F: data in frequencies
		yfiles = sorted(glob.glob('*z.tmp'))
		for k in range(0,len(yfiles)-1):
			expt = levmpy.Experiment(yfiles[k+1], outsuffix='OUTY%02d' % (k+1), path=os.getcwd())
			expt.fit()
			fp.loc[:,k] = expt.x
			fd.loc[:,k] = expt.freq
			yeee.loc[:,k] = expt.y1 + expt.y2*1j #This is the experimental admittance, it was automatically transformed by the LEVM library
		fp.to_csv('FitParametersZYRRF.csv', index=False) #main fitting function

def f(n, rs, rct, cdl, cad, alpha):
   return (rs + (1/((((n*1j)**alpha)*cdl)+(1/(rct+(1/(n*1j*cad))))))) #impedance fast plot function

def plot_fits():
	options = {'z':'z', 'Z':'Z', 'y':'y', 'Y':'Y'}
	selection = input("Choose Z for impedance plots or Y for admittance plots ")
	ret = options.get(selection)
	rs = fp.iloc[0,:]
	rct = fp.iloc[1,:]
	cdl = fp.iloc[2,:]
	cad = fp.iloc[3,:]
	alpha = fp.iloc[4,:]
	SMALL_SIZE = 12
	MEDIUM_SIZE = 14
	BIGGER_SIZE = 16
	plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
	plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
	plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
	plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
	plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
	plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
	plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
	plt.rc('figure', figsize=(8,8))
	if not ret:
		return plot_fits() #restart if wrong mode given
	if ret == 'z' or ret == 'Z':
		print('Impedance')
		for m in range(0,len(fp.iloc[0,:])): #for every file
			x = np.linspace(0.1, fd.iloc[0,m], 100000) #model values
			zr = f(x, float(rs.iloc[m]), float(rct.iloc[m]), float(cdl.iloc[m]), float(cad.iloc[m]), float(alpha.iloc[m])).real #function from parameters
			zi = -f(x, float(rs.iloc[m]), float(rct.iloc[m]), float(cdl.iloc[m]), float(cad.iloc[m]), float(alpha.iloc[m])).imag
			plt.plot(zr/1000, zi/1000)
			plt.scatter(zexp.iloc[:,m].apply(np.real)/1000, -zexp.iloc[:,m].apply(np.imag)/1000) #experimental scatter plots
		plt.xlabel(r'$Z_{real}\,/\,k\Omega$')
		plt.ylabel(r'$-Z_{imag}\,/\,k\Omega$')
		ah = max(zexp.apply(np.real).max().max()/1000, -zexp.apply(np.imag).min().min()/1000) #highest value of impedance real or imaginary
		plt.axis([0, ah, 0, ah]) #square axis for Nyquist
		plt.savefig('Z_plot.svg') #vector based output
		options = {'n':'n', 'N':'N', 'y':'yes', 'Y':'yes'}
		selection = input("Show plot? (y/n) ") #Interactive plot query
		ret = options.get(selection)
		if ret == 'yes' or ret == 'yes':
			plt.show()
		else:
			return print('...Interactive plot skipped')
	if ret == 'y' or ret == 'Y':
		print('Admittance')
		for m in range(0,len(fp.iloc[0,:])):
			x = np.linspace(0.1, fd.iloc[0,m], 100000)
			yr = (1/f(x, float(rs.iloc[m]), float(rct.iloc[m]), float(cdl.iloc[m]), float(cad.iloc[m]), float(alpha.iloc[m]))).real
			yi = (1/f(x, float(rs.iloc[m]), float(rct.iloc[m]), float(cdl.iloc[m]), float(cad.iloc[m]), float(alpha.iloc[m]))).imag
			plt.plot(yr*1e3, yi*1e3)
			plt.scatter(yeee.iloc[:,m].apply(np.real)*1e3, yeee.iloc[:,m].apply(np.imag)*1e3)
		plt.xlabel(r'$Y_{real}\,/\,mS$')
		plt.ylabel(r'$-Y_{imag}\,/\,mS$')
		ah = max(yeee.apply(np.real).max().max()*1e3, yeee.apply(np.imag).min().min()*1e3)
		plt.axis([0, ah, 0, ah])
		plt.savefig('Y_plot.svg')
		options = {'n':'n', 'N':'N', 'y':'yes', 'Y':'yes'}
		selection = input("Show plot? (y/n) ")
		ret = options.get(selection)
		if ret == 'yes' or ret == 'yes':
			plt.show()
		else:
			return print('...Interactive plot skipped') #main plotting function
fit_eqs()
plot_fits()


rfiles = glob.glob('*z.tmp')
for i in rfiles:
	os.remove(i)    #removes the ZYRRF files after use
