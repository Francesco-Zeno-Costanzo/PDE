import numpy as np
import matplotlib as mp
import scipy.sparse as sp
import scipy.sparse.linalg
import matplotlib.pyplot as plt
import matplotlib.animation as animation

N = 100
x = np.linspace(0, 1, N)
tstep = 5000
T = np.zeros((N, tstep))

#Temperatura iniziale
T_0 = 500*np.exp(-((x-0.5)/0.1)**2)

D = 0.2
dx = 0.01
dt = 1e-4
r = D*dt/dx**2
print(r)

A = sp.diags([-r, 1+2*r, -r], [-1, 0, 1], shape=(N, N)).toarray()

T[:, 0] = T_0
b = T_0

for i in range(tstep):

    u = np.linalg.solve(A, b)
    b = u
    T[:, i] = b

fig = plt.figure(1)
ax = fig.add_subplot(projection='3d')
gridx, gridy = np.meshgrid(range(tstep), x)
ax.plot_surface(gridx, gridy, T, cmap=mp.cm.coolwarm,vmax=250,linewidth=0,rstride=2, cstride=100)
ax.set_title('Diffusione del calore')
ax.set_xlabel('Tempo')
ax.set_ylabel('Lunghezza')
ax.set_zlabel('Temperatura')

fig = plt.figure(2)
plt.xlim(np.min(x), np.max(x))
plt.ylim(np.min(T), np.max(T))

line, = plt.plot([], [], 'b')
def animate(i):
    line.set_data(x, T[:,i])
    return line,


anim = animation.FuncAnimation(fig, animate, frames=tstep, interval=10, blit=True, repeat=True)

plt.grid()
plt.title('Diffusione del calore')
plt.xlabel('Distanza')
plt.ylabel('Temperatura')

#anim.save('calore.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

plt.show()