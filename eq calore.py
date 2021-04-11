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

tau=10
D=0.5*tau
dx=0.01
dt=1e-5 #dt<=dx**2/D
r = D*dt/dx**2

for time in range(1,tstep):
    for i in range(1,N-1):
        T[i,time]=T[i,time-1] + r*(T[i-1,time-1]+T[i+1,time-1]-2*T[i,time-1])
#T[0,time]=T[1,time] #per avere bordi non fissi
#T[N-1,time]=T[N-2,time]

fig = plt.figure()
ax = fig.gca(projection='3d')
gridx , gridy = np.meshgrid(range(tstep), range(N))
ax.plot_surface(gridx,gridy,T, cmap=mp.cm.coolwarm,vmax=250,linewidth=0,rstride=2, cstride=100)
ax.set_title('Diffusione del calore')
ax.set_xlabel('Tempo')
ax.set_ylabel('Lunghezza')
ax.set_zlabel('Temperatura')
plt.show()
