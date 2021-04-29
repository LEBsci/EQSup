import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import interpolate
from numpy.polynomial import Polynomial as p
'''
This test will be done based on potentials, maybe it would be a better choice
to do it based on current sign. Maybe not. Test is test.
1st step: text data is now ordered by potential from minimun to maximum.
'''
cycle_column = 10 # Selection of column, may change
area = 0.03952
data = pd.read_csv('Pt232222.csv', header=None, skiprows=1) # Read Ivium type data
cvraw = data.iloc[:,cycle_column:cycle_column+2].dropna() # Selects the working E, I columns and removes NaN values
emin = cvraw[cycle_column].idxmin() # Find the index for the minimun value
cv = pd.concat([cvraw.iloc[emin:,:],cvraw.iloc[0:emin,:]]).reset_index(drop=True) # Concatenation in order to later separate forward and reverse sweeps
esup = cv[cycle_column].idxmax()
cvpos = cv.iloc[1:esup,:].reset_index(drop=True) # Ascending order and fixing indexes for operations
cvneg = cv.iloc[esup+1:,:].sort_values(by=[cycle_column]).reset_index(drop=True) # Descending order and same as before
cvneg.iloc[:,1] = -cvneg.iloc[:,1] # Change current to positive values to later average
'''
2nd step: Obtain average cyclic voltammetry
'''
cvm = (cvpos + cvneg)/2
plt.scatter(cvm.iloc[:,0], cvpos.iloc[:,1]/area) # Two plots to select the adequate potential of hydrogen adsorption at terraces
plt.scatter(cvm.iloc[:,0], cvneg.iloc[:,1]/area) # An intersection point must be selected from which hydrogen adsorbs into de surface
plt.xlabel('$E$/V vs Pd-H')
plt.ylabel(r'$j\,/\,\mathrm{\mu A\,cm^{-2}}$')
plt.show()

ehmin = 0.268 # float(input('Write the starting potential of hydrogen adsorption (up to three decimals): '))
eint = np.linspace(cvm.iloc[0,0], ehmin, 1000)
ehmin_index = cvm.iloc[:,0].index[cvm.iloc[:,0] == ehmin].tolist()

dl1i = cvm.iloc[:,0].index[cvm.iloc[:,0] == 0.400].tolist() # index of 0.400 for double layer average
dl2i = cvm.iloc[:,0].index[cvm.iloc[:,0] == 0.499].tolist()
idl = np.average(cvm.iloc[dl1i[0]:dl2i[0]+1,1]) # average double layer current
print(idl/area)
plt.scatter(cvm.iloc[0:ehmin_index[0],0], (cvm.iloc[0:ehmin_index[0],1]-idl)/area) # Plot showing the average hydrogen adsorption at terraces minus the double layer current
plt.xlabel('$E$/V vs Pd-H')
plt.ylabel(r'$j\,/\,\mathrm{\mu A\,cm^{-2}}$')
plt.show()
'''
3rd is to integrate the charge via spline interpolation
'''
x = cvm.iloc[0:ehmin_index[0],0]
y = (cvm.iloc[0:ehmin_index[0],1]-idl)/area
xnew = np.linspace(cvm.iloc[0,0], cvm.iloc[ehmin_index[0]-1,0], 1000) # Defines the potentials for interpolation
tck = interpolate.splrep(x, y, s=0) # interpolated object of the function
ynew = interpolate.splev(xnew, tck, der=0) # interpolation evaluation if needed
qH2 = pd.DataFrame()
nu = 0.050
for n in range(0,len(xnew)):
    qH2.loc[n,0] = (interpolate.splint(cvm.iloc[0,0], xnew[n]+0.001, tck))
plt.plot(x, y, 'o', xnew, ynew, '-')
plt.xlabel('$E$/V vs Pd-H')
plt.ylabel(r'$j\,/\,\mathrm{\mu A\,cm^{-2}}$')
plt.show()
'''
4th step is to define the theoretical curve via Frumkin isotherm and compare
'''
T = 273.15
F = 96485
R = 8.31
e = 1.60217646e-13
d = 2.77e-8
n = 45 # base atoms in Pt(23,22,22), can later be adapted
s = (math.sqrt(3)/2)*(d**2)*(n-(1/3))*(math.sqrt(9*n**2+6*n+9)/(3*n-1)) # Inclined surface for electrode
qML = (n-1)*e/s-60
theta = qH2/(nu*qML) # hydrogen coverage coefficient
EH = pd.DataFrame(xnew).iloc[::-1].reset_index(drop=True) # converts xnew to dataframe for compatibility
dG0 = -F*EH-R*T*np.log(theta/(1-theta)) # Standard gibbs energy for compatibility
tM = theta.iloc[:,0].sub(0.35).abs().idxmin() # find values closest to desired number
tm = theta.iloc[:,0].sub(0.25).abs().idxmin() # this is done to obtain the linear regression and finding lim theta->0 dG0
z2 = p.fit(theta.iloc[tm:tM,0], dG0.iloc[tm:tM,0], 1)
plt.scatter(theta.iloc[:,0],dG0.iloc[:,0])
plt.scatter(theta.iloc[tm:tM,0],dG0.iloc[tm:tM,0])
plt.plot(np.linspace(0,0.5,1000),z2(np.linspace(0,0.5,1000)), c='black')
plt.show()
r = z2.deriv()
g = r(0)/(R*T)
dG0i = z2(0)
thetat = np.linspace(0.001, 0.999, 999)
Et = -dG0i/F-(R*T/F)*np.log(thetat/(1-thetat))-r(0)*thetat/F
jt = qML*nu*F/(R*T)*((thetat*(1-thetat))/(1+thetat*g*(1-thetat)))
plt.scatter(cvm.iloc[0:ehmin_index[0],0], (cvm.iloc[0:ehmin_index[0],1]-idl)/area)
plt.plot(Et,jt)
plt.show()
print(r(0), g, qML)
