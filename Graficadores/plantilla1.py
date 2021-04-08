import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))    #Detecta la ubicacion del archivo
data = os.path.join(THIS_FOLDER, 'bigdata.csv') #Consigue el archivo de datos para graficar. Ver archivo para comparar el formato

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


pt = pd.read_csv(data, index_col='electrode')   #La columna 'electrode' se usa para identificar. Puede tener otra informacion
miller = pt.index.array.unique().dropna()   #Esto hace una lista de cuantos 'electrodos' hay

for n in miller:    #Ciclo usando los valores de miller
    plt.plot(pt.loc[n].iloc[:,0], pt.loc[n].iloc[:,5], label=n) #En el ciclo se grafica cada uno de los 'electrodos' y se guardan
plt.axhline(0, 0, 1, c = 'blue')   #Se grafica el eje cero
plt.xlabel('$E$/V vs Pd/H') #Nombre de eje X
plt.ylabel(r'$i\,/\,\mathrm{\mu A}$')   #Nombre de eje Y
plt.legend()    #Se ubica la leyenda en el mejor lugar automaticamente
plt.savefig('voltagramas.svg')  #Se guarda la grafica en formato de vectores
plt.show()  #Se muestra el grafico de forma interactiva
