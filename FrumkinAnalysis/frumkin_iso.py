import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.integrate as spi

'''
This test will be done based on potentials, maybe it would be a better choice
to do it based on current sign. Maybe not. Test is test.
1st step: text data is now ordered by potential from minimun to maximum.
'''

cycle_column = 8 # Selection of column, may change
data = pd.read_csv('Pt232222.csv', header=None, skiprows=1) # Read Ivium type data
cvraw = data.iloc[:,cycle_column:cycle_column+2].dropna() # Selects the working E, I columns and removes NaN values
emin = cvraw[cycle_column].idxmin() # Find the index for the minimun value
cv = pd.concat([cvraw.iloc[emin:,:],cvraw.iloc[0:emin,:]]).reset_index(drop=True) # Concatenation in order to later separate forward and reverse sweeps
esup = cv[cycle_column].idxmax()
cvpos = cv.iloc[1:esup,:].reset_index(drop=True) # Ascending order and fixing indexes for operations
cvneg = cv.iloc[esup+1:,:].sort_values(by=[cycle_column]).reset_index(drop=True) # Descending order and same as before
cvneg.iloc[:,1] = -cvneg.iloc[:,1] # Change current to positive values to later average

'''
2nd step: Interpolate and obtain average cyclic voltammetry
'''
cvm = (cvpos + cvneg)/2
plt.scatter(cvm.iloc[:,0], cvpos.iloc[:,1]) # Two plots to select the adequate potential of hydrogen adsorption at terraces
plt.scatter(cvm.iloc[:,0], cvneg.iloc[:,1])
plt.show()

ehmin = float(input('Write the starting potential of hydrogen adsorption (up to three decimals): '))
# eint = np.linspace(cvm.iloc[0,0], ehmin, 1000)
ehmin_index = cvm.iloc[:,0].index[cvm.iloc[:,0] == ehmin].tolist()
# cvHt = cvm.iloc[:ehmin,:]
# i = np.arange(y, cvm.iloc[0:,0])
plt.scatter(cvm.iloc[0:ehmin_index[0],0], cvm.iloc[0:ehmin_index[0],1]) # Plot showing the average hydrogen adsorption at terraces
plt.show()
# print(ehmin_index[0])
# cvraw.to_csv('out.csv')
