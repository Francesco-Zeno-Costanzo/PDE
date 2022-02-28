import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

N = 100
t = 500
c = 0.3
y = np.linspace(0, 1, N)

U = np.zeros((N, N, t))
U[:,:,0] = np.exp(-((y - 0.5)**2)/0.05)
U[:,:,0] = -U[:,:,0]*np.transpose(U[:,:,0])

for j in np.arange(0,t-1):
    for i in np.arange(0,N-1):
        for k in np.arange(0,N-1):
            U[i,k,j+1] = 2*U[i,k,j]-U[i,k,j-1]+ c*(U[i+1,k,j]-2*U[i,k,j]+U[i-1,k,j]+U[i,k+1,j]-2*U[i,k,j]+U[i,k-1,j])

##
#la simulazione restituisce un'onda la cui ampiezza e maggiore della
#condizione inizale. Probabilmente ciò è dovuto a problemi del metodo adottato
##

##grafico
h = 130 #istante del tempo a cui vedere il grafico
ZM = np.max(U)
zm = np.min(U)

fig = plt.figure(1)
ax = fig.gca(projection='3d')

ax.set_xlabel('Distanza')
ax.set_ylabel('Distanza')
ax.set_zlabel('Ampiezza')


gridx, gridy = np.meshgrid(y, y)
ax.plot_surface(gridx,gridy,U[:,:,h], cmap=cm.coolwarm)

plt.show()