import numpy as np
import matplotlib as mp
import matplotlib.pyplot as plt
import matplotlib.animation as animation

par = np.loadtxt(r'input_c.txt')
N = int(par[0])
t = int(par[1]) + 1
D = par[2]
dt = par[3]
dx = par[4]

u = np.loadtxt(r'calore.dat')

T = np.reshape(u, (t, N))

x = np.linspace(0, N*dx, N)


fig = plt.figure(1)
ax = fig.add_subplot(projection='3d')
gridx, gridy = np.meshgrid(x, range(t))
ax.plot_surface(gridx, gridy, T, cmap=mp.cm.coolwarm,vmax=25,linewidth=0,rstride=2, cstride=2)
ax.set_title('Diffusione del calore')
ax.set_xlabel('Tempo')
ax.set_ylabel('Lunghezza')
ax.set_zlabel('Temperatura')

fig = plt.figure(2)
plt.xlim(np.min(x), np.max(x))
plt.ylim(np.min(T), np.max(T))

line, = plt.plot([], [], 'b')
def animate(i):
    line.set_data(x, T[i,:])
    return line,


anim = animation.FuncAnimation(fig, animate, frames=t, interval=10, blit=True, repeat=True)

plt.grid()
plt.title('Diffusione del calore')
plt.xlabel('Distanza')
plt.ylabel('Temperatura')

plt.show()
