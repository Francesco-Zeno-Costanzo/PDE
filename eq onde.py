import numpy as np
import matplotlib as mp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation


N = 125
t = 315
c = 1
y = np.linspace(0, 1, N)
U = np.zeros((N,t))
U[0:N,0] = np.sin(2*np.pi*y)

for j in np.arange(0,t-1):
    for i in np.arange(0,N-1):
        U[i, j+1] = 2*U[i,j]-U[i,j-1]+ c*(U[i+1, j] - 2*U[i, j] + U[i-1, j])

##
#la simulazione restituisce un'onda la cui ampiezza e maggiore della
#condizione inizale. Probabilmente ciò è dovuto a problemi del metodo adottato
##

fig = plt.figure(1)
ax = fig.gca(projection='3d')

ax.set_title('Propagazione onde')
ax.set_xlabel('Tempo')
ax.set_ylabel('Lunghezza')
ax.set_zlabel('Ampiezza')


gridx, gridy = np.meshgrid(range(t), y)
ax.plot_surface(gridx,gridy,U)

fig = plt.figure(2)
plt.xlim(np.min(y), np.max(y))
plt.ylim(np.min(U)-0.1, np.max(U)+0.1)

line, = plt.plot([], [], 'b')
def animate(i):
    line.set_data(y, U[:,i])
    return line,


anim = animation.FuncAnimation(fig, animate, frames=t ,interval=20, blit=True, repeat=True)

plt.grid()
plt.title('Propagazione onde')
plt.xlabel('Distanza')
plt.ylabel('Ampiezza')

#anim.save('onda.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

plt.show()