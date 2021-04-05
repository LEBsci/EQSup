import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

cy = 8 # Selection of column, may change
data = pd.read_csv('Pt322.csv', header=None, skiprows=1) # Read Ivium type data
cvraw = data.iloc[:,cy:cy+2].dropna() # Selects the working E, I columns and removes NaN values
emin = cvraw[cy].idxmin() # Find the index for the minimun value
cv = pd.concat([cvraw.iloc[emin:,:],cvraw.iloc[0:emin,:]]) # Concatenation in order to later separate forward and reverse sweeps
esup = cvraw[cy].idxmax()
cvpos = cv.loc[:esup,:]
cvneg = cv.loc[esup:,:]
print(cvpos)

# cvraw.to_csv('out.csv')
plt.scatter(cv.iloc[:,0], cv.iloc[:,1])
plt.show()
'''
version 0.1.2
This test will be done based on potentials, maybe it would be a better choice
to do it based on current sign. Maybe not. Test is test.
text data is now ordered as intended for one cycle
'''
