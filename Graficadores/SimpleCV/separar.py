"""Método para separar los datos de corriente positiva y negativa en la voltametría cíclica y posteriormente graficar
Basado principalmente en una Boolean mask para filtrar los datos de un array"""

import numpy as np
import matplotlib.pyplot as plt
import os

THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))    #Detecta la ubicacion del archivo
data = os.path.join(THIS_FOLDER, 'Pt554_25C.txt') #Consigue el archivo de datos para graficar. Ver archivo para comparar el formato

def plotfunction(cvs):
    plt.scatter(cvs[:,0], cvs[:,1]*10e6, c='r')
    plt.axhline(0, 0, 1, c = 'black')   #Se grafica el eje cero
    plt.xlabel('$E$/V vs Pd/H') #Nombre de eje X
    plt.ylabel(r'$i\,/\,\mathrm{\mu A}$')   #Nombre de eje Y
    plt.show()

cv = np.loadtxt(data, skiprows=14)
plt.title('Pt(554)')
# plt.savefig('Pt554.png', dpi=300)
plotfunction(cv)


cvNeg = cv[cv[:,1] < 0] #Boolean mask para encontrar los valores de corriente negativa
np.savetxt("negative.txt", cvNeg)
plt.title('Pt(554) Negative')
# plt.savefig('Pt554.png', dpi=300)
plotfunction(cvNeg)

cvPos = cv[cv[:,1] > 0]
np.savetxt(os.path.join(THIS_FOLDER, 'positive.txt'), cvPos)
plt.title('Pt(554) Positive')
# plt.savefig('Pt554.png', dpi=300)
plotfunction(cvPos)