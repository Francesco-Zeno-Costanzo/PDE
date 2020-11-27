import numpy as np
import matplotlib as mp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

N=125
t=15
c=1
y=np.linspace(0, N, N)
U=np.zeros((N,t))
U[:,0]=np.sin(0.05*y)
for j in np.arange(0,t-2):
    for i in np.arange(1,N):
        U[i,j] = c*(2*U[i,j]-U[i,j+2]-U[i,j-2])

fig = plt.figure(1)
ax = fig.gca(projection='3d')
gridx , gridy = np.meshgrid(range(t),range(N))
ax.plot_surface(gridx,gridy,U)
ax.set_title('Propagazione onde')
ax.set_xlabel('Tempo')
ax.set_ylabel('Lunghezza')
ax.set_zlabel('Ampiezza')
plt.show()