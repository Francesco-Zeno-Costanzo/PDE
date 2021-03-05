import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

N=100
t=500
c=0.5
y=np.linspace(0, N, N)

U=np.zeros((N,N,t))
U[:,:,0]=np.exp(-((y-N//2)**2)/100)
U[:,:,0]=-U[:,:,0]*np.transpose(U[:,:,0])

for j in np.arange(0,t-1):
    for i in np.arange(0,N-1):
        for k in np.arange(0,N-1):
            U[i,k,j+1] = 2*U[i,k,j]-U[i,k,j-1]+ c*(U[i+1,k,j]-2*U[i,k,j]+U[i-1,k,j]+U[i,k+1,j]-2*U[i,k,j]+U[i,k-1,j])

##grafico
h=10 #istante del tempo a cui vedere il grafico
ZM=np.max(U)
zm=np.min(U)

fig = plt.figure(1)
ax = fig.gca(projection='3d')

ax.set_xlabel('Distanza')
ax.set_ylabel('Distanza')
ax.set_zlabel('Ampiezza')
ax.set_zlim(zm, ZM)

gridx , gridy = np.meshgrid(range(N), range(N))
ax.plot_surface(gridx,gridy,U[:,:,h])

plt.show()