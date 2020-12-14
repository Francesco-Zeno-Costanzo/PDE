import numpy as np
import matplotlib as mp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

N=125
t=314
c=1
y=np.linspace(0, N, N)
U=np.zeros((N,t))
U[0:N,0]=np.sin(0.05*y)

for j in np.arange(0,t-2):
    for i in np.arange(0,N-1):
        U[i, j+1] = 2*U[i,j]-U[i,j-1]+ c*(U[i+1, j]-2*U[i, j]+U[i-1, j])

fig = plt.figure(1)
ax = fig.gca(projection='3d')

ax.set_title('Propagazione onde')
ax.set_xlabel('Tempo')
ax.set_ylabel('Lunghezza')
ax.set_zlabel('Ampiezza')


gridx , gridy = np.meshgrid(range(t), range(N))
ax.plot_surface(gridx,gridy,U)

plt.show()