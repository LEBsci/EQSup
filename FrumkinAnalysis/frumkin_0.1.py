import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

cy = 8
data = pd.read_csv('Pt322.csv', header=None, skiprows=1)
cvraw = data.iloc[:,cy:cy+2].dropna()
cv = cvraw.sort_values(by=[cy])
# plt.scatter(cv.iloc[:,0], cv.iloc[:,1])
# plt.show()
print(cv)
# cv.to_csv('out.csv')
'''
version 0.1
ready to plot one cv as a scatter plot but it is not ordered
'''
