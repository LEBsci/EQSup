import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import splrep, splev, splint
import glob, os, random
from matplotlib import cm
from matplotlib.colors import Normalize

'''Plot settings'''

SMALL_SIZE = 16
MEDIUM_SIZE = 18
BIGGER_SIZE = 20

#plt.style.use('bmh')

plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=SMALL_SIZE)  # fontsize of the figure title
plt.rc('figure', figsize=(8,8))

def integ(x, i, tck, constant=-1):
    out = np.zeros(x.shape)
    sup = x[i == 0]
    top = np.where(x == sup[0])[0][0]
    for n in range(top):
        out[n] = splint(x[0], x[n], tck)
    # out += constant
    return out

def isogQ(cycle):
    ifiles = sorted(glob.glob("*C.txt")) #First step import data files
    q = np.zeros(len(ifiles))
    qx = np.zeros((len(ifiles),10000))

    for i in range(len(ifiles)): #will do a for cycle for the whole system
        data=np.loadtxt(ifiles[i], skiprows=1)
        

        cvpos = data[(data[:,2*cycle-1] > 0)] #Select positive currents 
        limpos = np.argmax(cvpos[:,2*(cycle-1)]) #Find maximum potential to discard duplicates
        

        cvneg = data[(data[:,2*cycle-1] < 0)]
        limneg  = np.argmin(cvneg[:,2*(cycle-1)]) #Find minimum potential to discard duplicates
        limnegp = np.argmax(cvneg[:,2*(cycle-1)])
        eneg = np.concatenate((cvneg[limneg::-1,2*(cycle-1)],cvneg[-1:limnegp:-1,2*(cycle-1)])) #organize lower values since the scan started at 0.1
        ineg = np.concatenate((cvneg[limneg::-1,2*cycle-1],cvneg[-1:limnegp:-1,2*cycle-1]))

        tckp = splrep(cvpos[:limpos,2*(cycle-1)], cvpos[:limpos,2*cycle-1], s=0) #Positive interpolation
        tckn = splrep(eneg, ineg, s=0)
    
        xpos = np.linspace(cvpos[0,2*(cycle-1)], cvpos[limpos,2*(cycle-1)], 10000) #new potential values for interpolation
        xneg = np.linspace(eneg[0], eneg[-1], 10000)
        xavg = np.linspace(max(xpos[0], xneg[0]), min(xpos[-1], xneg[-1]), 10000) #potential values for average results, selecting the appropriate intervale in which both exist

        
        
        ipos = splev(xpos, tckp, der=0)
        ineg = splev(xneg, tckn, der=0)
        iavg = (ipos-ineg)/2
        inew = iavg-min(iavg[(xavg>0.3)*(xavg<0.4)])

        tcki = splrep(xavg, inew, s=0)
        q[i] = splint(xavg[0], xavg[inew == 0], tcki)/0.05 #integration from the lowest value to the start of hydrogen adsorption
        qx[i] = integ(xavg, inew, tcki)/0.05
        plt.plot(xavg, qx[i])
    plt.show()      
    return q, qx

qHNP = 230
qH, qe = isogQ(3)
areas = qH/(230)
print(qe[0])

# def isogt():
#     ifiles = sorted(glob.glob("*C.txt")) #First step import data files
#     for i in range(len(ifiles)):
#         thetaH[i] = qH[i] 

'''


area=0.05227
F=96485
R=8.31
qML=180 # Esto debería ser 240, pero ajusta mejor usando 170
v=0.05
data=np.loadtxt('voltametriasm.dat',skiprows=2)
# en filas tengo distintos potenciales y en columnas distintas temperaturas
E=data[:,0] # La primera columna es la de potenciales  
limsup=np.argmax(E)   #Busco el indice correspondiente al limite superior 
Epos=E[:limsup]  #Separo barrido positivo y negativo
Eneg=E[limsup:]
Eint=np.linspace(0.01,0.79,1000)  #Esto lo usaré para interpolar
EHmin=0.48 #Este es el potencial al que empieza a adsorberse H
EH=Eint[(Eint<EHmin)] #Extraigo los potenciales para las curvas de H
imaxH=len(EH)-1
EOHmin=0.49  #Este es el potencial al que empieza a adsorberse OH
EOH=Eint[(Eint>EOHmin)] #Extraigo los potenciales para las curvas de OH
iminOH=np.argwhere(Eint>EOHmin)[0,0]
corr=data[:,1:] #extraigo las corrientes
temp=np.array([10,15,20,25,30,35,40])
temp+=273  
cols=len(temp)
corrint=np.zeros((1000,cols))  #inicializo una matriz para contener las corrientes interpoladas
qH2=np.zeros((len(EH),cols))  #inicializo qH2
qOH2=np.zeros((len(EOH),cols))  #inicializo qOH2

cmap = cm.get_cmap('plasma', cols)
norm = Normalize(vmin=np.min(0),vmax=np.max(cols)) 

for i in range(0,int(corr.shape[1]/3)):
    j=corr[:,i*3+1]-0.02 #hay tres ciclos, me quedo con el segundo
    plt.plot(E,j/area,label=temp[i],color=cmap(norm(i)))
    jpos=j[0:limsup]
    jneg=-j[limsup:]
    jposint=interp1d(Epos,jpos)(Eint)  #interpolo las corrientes para promediar
    jnegint=interp1d(Eneg,jneg)(Eint)
    jmed=(jposint+jnegint)/2  #promedio
    jdl=jmed[(Eint>0.4)*(Eint<0.6)]   #busco el valor de doble capa   
    jcorr=jmed-min(jdl)  # corriente media corregida la doble capa   
    corrint[:,i]=jcorr   #almaceno el resultado en corrint
    q=cumtrapz(jcorr,Eint,initial=0)/v    #integro para hallar la carga
    qH=q[imaxH]-q[:imaxH+1]     #separo la curva de hidrógeno
    qH2[:,i]=qH      #y la almaceno en qH2
    qOH=q[iminOH:]-q[iminOH]    #separo la curva de OH
    qOH2[:,i]=qOH   #y la almaceno en qOH2
 

plt.axhline(0, 0, 1, c = 'black')
plt.xlabel('$E$/V vs RHE')
plt.ylabel(r'$j\,/\,\mathrm{\mu A \,cm^{-2}}$')
plt.legend(title='Temperatures')
plt.savefig('voltagramas_sosa.svg')
# plt.show()
plt.close()
for i in range(0,cols):
    plt.plot(Eint,corrint[:,i]/area,label=temp[i])
plt.xlabel('$E$/V vs RHE')
plt.ylabel(r'$j\,/\,\mathrm{\mu A \,cm^{-2}}$')
plt.legend(title='temperature / K')
plt.savefig('voltagramaspromedio_sosa.svg')
# plt.show()    
plt.close()
#area=np.max(qH2)/140
thetaH=qH2/area/qML
charge = pd.DataFrame((qH2[1,:],qOH2[-1,:]))
charge.to_csv('carga.csv', index=False)

for i in range(0,cols):
    plt.plot(EH,qH2[:,i]/area,label=temp[i])

plt.legend(title='temperature / K')
plt.xlabel('$E$/V vs RHE')
plt.ylabel(r'$q\,/\,\mathrm{\mu C \,cm^{-2}}$')
plt.xlim(0,0.5)
plt.ylim(0,np.max(qH2)*1.1/area)
plt.twinx()

for i in range(0,cols):
    plt.plot(EH,thetaH[:,i])
plt.ylabel(r'$\theta$')    
plt.ylim(0,np.max(qH2)*1.1/area/qML)

plt.savefig('cargasH_sosa.svg')
# plt.show()
plt.close()
thetaint=np.linspace(0.0001,0.9999,1000)
Q=np.log(thetaH/(1-thetaH))
DeltaG=np.zeros_like(Q)
r=[]
DG0=[]               #inicializo

for i in range(0,cols):
    DeltaG[:,i]=-F*EH-R*temp[i]*Q[:,i]   
    filtro=(thetaH[:,i]>0.1)*(thetaH[:,i]<0.5) # busco los calores entre 0.1 y  0.5
    DG=DeltaG[filtro,i]
    tH=thetaH[filtro,i]
    z=np.polyfit(tH,DG,1)  #ajuste lineal
    ri,DG0i=z      # esto me sirve para simplificar la notación después
    DG0.append(DG0i)  # ordenada en el origen 
    r.append(ri)   # pendiente 
    #calculo la curva teorica
    Ec=-DG0i/F-R*temp[i]/F*np.log(thetaint/(1-thetaint))-ri/F*thetaint
    corrcalc=F/R/temp[i]*qML*v*\
    thetaint*(1-thetaint)/(1+ri/R/temp[i]*thetaint*(1-thetaint))
    if i==0 or i==cols-1:
        plt.plot(Ec,corrcalc,'--',color='C%i'%i) 
        plt.plot(Eint,corrint[:,i]/area,'-',color='C%i'%i, label=temp[i]) #solo represento la primera y la última
plt.xlabel('$E$/V vs RHE')
plt.ylabel(r'$j\,/\,\mathrm{\mu A \,cm^{-2}}$')
plt.legend(title='temperature / K')
plt.xlim(-0.10,0.5)
plt.ylim(0,40)
plt.savefig('corrientesH_calc_sosa.svg')
# plt.show()
plt.close()

 
for i in range(0,cols):
    plt.plot(thetaH[:,i],-DeltaG[:,i]/1000,label=temp[i], color=cmap(norm(i)))
plt.xlabel(r'$\theta_{H}$')
plt.ylabel(r'$-\Delta G_{ads}\,/\,\mathrm{kJ\, mol^{-1}}$') 
plt.legend(title='temperature / K')

dx = np.zeros_like(DeltaG)
dy = np.zeros_like(thetaH)
d = dx/dy
for i in range(0,cols): #calculo de derivadas para obtener parametro de interaccion
    dy[:,i] = np.gradient(-DeltaG[:,i]/1000)
    dx[:,i] = np.gradient(thetaH[:,i])
    d[:,i] = dy[:,i]/dx[:,i]
param = pd.DataFrame(d)
param.to_csv('frumkinpara2.csv', index=False)
theta = pd.DataFrame(thetaH)
theta.to_csv('frumkintheta2.csv', index=False)
plt.savefig('DeltaGH_sosa.svg')
# plt.show()
plt.close()

# 
# 
# plt.plot(temp,r,'o',label='r/J/mol')
# plt.ylim(0,23e3)     
# plt.legend()
# plt.xlabel(r'$T$/K')
# plt.ylabel(r'$r/\mathrm{J \,mol^{-1}}$') 
# 
# ax.twinx()
# plt.plot(temp,DG0,'^',label=r'$\Delta G^0/\mathrm{J \,mol^{-1}}$')
# plt.legend()
# plt.ylabel(r'$\Delta G^0/\mathrm{J \,mol^{-1}}$') 


#primero hago el cálculo con 10 recubrimientos para hacer la figura
thetaint=np.linspace(0.01,0.5,10)
filas=len(thetaint)
DeltaGint=np.zeros((filas,cols))    #inicializo para guardar DeltaG interpolado
indice_ajuste=1      #   1 sería ajuste lineal, 2 a una parabola
for i in range(0,cols):
    DeltaGint[:,i]=interp1d(thetaH[:,i],DeltaG[:,i])(thetaint)
for i in range(0,filas):
    plt.plot(temp,-DeltaGint[i,:]/1000,'o',label=round(thetaint[i],2))
    z=np.polyfit(temp,DeltaGint[i,:],indice_ajuste)
    ffit=np.poly1d(z)
    xnew=np.linspace(min(temp),max(temp),100)
    ynew=ffit(xnew)
    plt.plot(xnew,-ynew/1000,'-') 
plt.xlabel(r'$T$/K')
plt.ylabel(r'$-\Delta G_{ads}\,/\,\mathrm{kJ\, mol^{-1}}$')
plt.legend(title='coverage')
plt.savefig('DeltaGH_vs_temp_sosa.svg')     
# plt.show()
plt.close()

#repito ahora con más puntos
thetaint=np.linspace(0.01,0.5,100)
filas=len(thetaint)
DeltaGint=np.zeros((filas,cols))    
for i in range(0,cols):
    DeltaGint[:,i]=interp1d(thetaH[:,i],DeltaG[:,i])(thetaint)
DeltaS=np.zeros((filas,cols))
for i in range(0,filas):
    z=np.polyfit(temp,DeltaGint[i,:],indice_ajuste)
#    DeltaS[i]=-z[0]
    ffit=np.poly1d(z)
    derivada=np.polyder(ffit,1)(temp)
    DeltaS[i]=-derivada
for i in range(0,cols):
    plt.plot(thetaint,DeltaS[:,i],label=temp[i], color='k')
entr = pd.DataFrame((thetaint,DeltaS[:,1])).T
entr.to_csv('dS.csv', index=False)
# plt.savefig('DeltaGH_sosa.svg')
plt.xlabel(r'$\theta_H$')
plt.ylabel(r'$\Delta S_{ads}\,/\,\mathrm{J \,mol^{-1}}$')

plt.savefig('DeltaSH_sosa.svg')  
# plt.show()
plt.close()


# np.savetxt('entropia_sosa.dat',np.vstack((thetaint,DeltaS)).transpose(),)
# np.savetxt('frumkin_sosa.dat',np.vstack((temp,DG0,r)).transpose(),)





DeltaH=np.zeros_like(DeltaGint)
for i in range (0,cols):
    DeltaH[:,i]=DeltaGint[:,i]+temp[i]*DeltaS[:,i]
    plt.plot(thetaint,DeltaH[:,i]/1000,label=temp[i], color=cmap(norm(i)))

plt.xlabel(r'$\theta_{H}$')
plt.ylabel(r'$\Delta H_{ads}\,/\,\mathrm{kJ\, mol^{-1}}$')
plt.legend(title='temperature/K')    

plt.savefig('DeltaHH_sosa.svg')
plt.close()
# plt.show()
DH2=436
EPtH=0.5*DH2-DeltaH
for i in range (0,cols):
    plt.plot(thetaint,EPtH[:,i]/100,label=temp[i], color=cmap(norm(i)))

plt.xlabel(r'$\theta_H$')
plt.ylabel(r'$E_{Pt-H}\,/\,\mathrm{kJ\, mol^{-1}}$')
plt.legend(title='temperature/K')    

plt.savefig('EPtH_sosa.svg')

# plt.show() 
plt.close()


plt.rc('font', size=28)          # controls default text sizes
plt.rc('axes', titlesize=28)     # fontsize of the axes title
plt.rc('axes', labelsize=28)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=28)    # fontsize of the tick labels
plt.rc('ytick', labelsize=28)    # fontsize of the tick labels
plt.rc('legend', fontsize=28)    # legend fontsize
plt.rc('figure', titlesize=28)  # fontsize of the figure title

dHT1=[DeltaH[-1,0]/-1000]
index=int(len(DeltaH[:,0])*(4/5))
dHT2=[DeltaH[index,0]/-1000]


for i in range(1,cols):
    dHT1.append(DeltaH[-1,i]/-1000)
    dHT2.append(DeltaH[index,i]/-1000)

plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = False
plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = True

fig, ax = plt.subplots()

plt.scatter(temp, dHT2, label=r'$\theta=$'+str(round(thetaint[index],1)), color='red')
plt.scatter(temp, dHT1, label=r'$\theta=0.5$', color='blue')
plt.legend()    
ax.set_title(r'$T$/K')
plt.ylabel(r'$-\Delta H_{ads}\,/\,\mathrm{kJ\, mol^{-1}}$')
plt.savefig('insert1.svg')
plt.show()

Ep1=[EPtH[-1,0]/100]
index=int(len(EPtH[:,0])*(4/5))
Ep2=[EPtH[index,0]/100]


for i in range(1,cols):
    Ep1.append(EPtH[-1,i]/100)
    Ep2.append(EPtH[index,i]/100)

plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = True
plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = False

plt.scatter(temp, Ep2, label=r'$\theta=$'+str(round(thetaint[index],1)), color='red')
plt.scatter(temp, Ep1, label=r'$\theta=0.5$', color='blue')
plt.legend()    
plt.xlabel(r'$T$/K')
plt.ylabel(r'$E_{Pt-H}\,/\,\mathrm{kJ\, mol^{-1}}$')
plt.savefig('insert2.svg')
plt.show()

#np.savetxt('results.dat',np.vstack(qH2))
'''