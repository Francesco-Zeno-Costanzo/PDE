import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

par = np.loadtxt(r'input_onde.txt')
N = int(par[0]) + 1
t = int(par[1]) + 1
D = par[2]
dt = par[3]
dx = par[4]

u = np.loadtxt(r'onde_lw.dat')

T = np.reshape(u, (t, N))

x = np.linspace(0, N*dx, N)


fig = plt.figure(1)
ax = fig.add_subplot(projection='3d')
gridx, gridy = np.meshgrid(x, range(t))
ax.plot_surface(gridx, gridy, T)
ax.set_title('Corda elastica')
ax.set_xlabel('distanza')
ax.set_ylabel('tempo')
ax.set_zlabel('Ampiezza')

fig = plt.figure(2)
plt.xlim(np.min(x), np.max(x))
plt.ylim(np.min(T), np.max(T))

line, = plt.plot([], [], 'b')
def animate(i):
    line.set_data(x, T[i,:])
    return line,


anim = animation.FuncAnimation(fig, animate, frames=t, interval=10, blit=True, repeat=True)

plt.grid()
plt.title('Corda elastica')
plt.xlabel('Distanza')
plt.ylabel('Ampiezza')

plt.show()
