from operator import index
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

scanrate = 0.05

data = np.loadtxt('Pt215.txt', delimiter=',',skiprows=1)

def separate(texto, cycle):
    cvneg = texto[texto[:,cycle] < 0]
    cvpos = texto[texto[:,cycle] > 0]
    


'''
THIS INTERPOLATE VERSION IS NOW OBSOLETE
'''



corr = data[:,1:] # extraigo todas las corrientes
EHmin=0.48 # Este es el potencial al que empieza a adsorberse H
Enew = np.arange(0.001, 0.790, 0.001) # valores teóricos de potencial
EvH = Enew[(Enew<EHmin)]
jnew = np.zeros((len(Enew),len(corr[0,:])))
qH2=np.zeros((len(corr[0,:])))  #inicializo qH2

for i in range(0,len(corr[0,:])):
    jpos = corr[0:limsup,i]/area
    jneg = -corr[limsup:,i]/area
    tckp = interpolate.splrep(Epos, jpos, s=0) # interpolación
    tckn = interpolate.splrep(Eneg[::-1], jneg[::-1], s=0) # debo usar Eneg de menor a mayor
    jposnew = interpolate.splev(Enew, tckp, der=0) # Corriente interpolada
    jnegnew = interpolate.splev(Enew, tckn, der=0)
    jmed = (jposnew+jnegnew)/2 # promedio
    jdl=jmed[(Enew>0.4)*(Enew<0.6)] # busco el valor de doble capa  
    jc=jmed-min(jdl) # corriente media corregida la doble capa
    jnew[:,i] = jc # almaceno el resultado
    intp = interpolate.splint(0.001, 0.480, tckp) # integral de la corriente interpolada
    intn = interpolate.splint(0.001, 0.480, tckn)
    qH2[i] = (intp+intn)/(scanrate*2)

np.savetxt('qrev.csv', qH2,delimiter=',')
