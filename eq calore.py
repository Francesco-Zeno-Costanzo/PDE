import numpy as np
import matplotlib as mp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

N=100
x=np.linspace(0, N, N)
tstep=5000
T=np.zeros((N,tstep))
#Temperatura iniziale
T[0:N,0]=500*np.exp(-((50-x)/20)**2)

ts=0.1
tau=0.9*ts


for time in range(1,tstep):
    for i in range(1,N-1):
        T[i,time]=T[i,time-1]+0.5*tau/ts*(T[i-1,time-1]+T[i+1,time-1]-2*T[i,time-1])
    T[0,time]=T[1,time]
    T[N-1,time]=T[N-2,time]

fig = plt.figure()
ax = fig.gca(projection='3d')
gridx , gridy = np.meshgrid(range(tstep),range(N))
ax.plot_surface(gridx,gridy,T,cmap=mp.cm.coolwarm,vmax=250,linewidth=0,rstride=2, cstride=100)
ax.set_title('Diffusione del calore')
ax.set_xlabel('Tempo')
ax.set_ylabel('Lunghezza')
ax.set_zlabel('Temperatura')
plt.show()