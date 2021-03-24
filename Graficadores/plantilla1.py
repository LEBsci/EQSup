import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))
data = os.path.join(THIS_FOLDER, 'bigdata.csv')

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


pt = pd.read_csv(data, index_col='electrode')
miller = pt.index.array.unique().dropna()

for n in miller:
    plt.plot(pt.loc[n].iloc[:,0], pt.loc[n].iloc[:,5], label=n)
plt.axhline(0, 0, 1, c = 'black')
plt.xlabel('$E$/V vs Pd/H')
plt.ylabel(r'$i\,/\,\mathrm{\mu A}$')
plt.legend()
plt.savefig('voltagramas.svg')
plt.show()
