import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Parameters
par = np.loadtxt(r'input_wave.txt')
N   = int(par[0]) + 1
t   = int(par[1]) + 1
D   = par[2]
dt  = par[3]
dx  = par[4]

# load solutions
u = np.loadtxt(r'wave_lw.dat')
# reshape for plot
T = np.reshape(u, (t, N))
# x axis
x = np.linspace(0, N*dx, N)

#==========================================================
# Plot
#==========================================================

fig = plt.figure(1)
ax = fig.add_subplot(projection='3d')
gridx, gridy = np.meshgrid(x, np.linspace(0, t*dt, t))
ax.plot_surface(gridx, gridy, T)
ax.set_title('Wave propagration')
ax.set_xlabel('Distance')
ax.set_ylabel('time')
ax.set_zlabel('Amplitude')

fig = plt.figure(3)
ax = fig.add_subplot(projection='3d')
gridx, gridy = np.meshgrid(x, np.linspace(0, t*dt, t))
sol = np.cos(2*np.pi*0.5*gridy)*np.sin(2*np.pi*gridx)
ax.plot_surface(gridx, gridy, sol)
ax.set_title('Wave propagration')
ax.set_xlabel('Distance')
ax.set_ylabel('time')
ax.set_zlabel('Amplitude')

fig = plt.figure(2)
plt.xlim(np.min(x), np.max(x))
plt.ylim(np.min(T), np.max(T))

line, = plt.plot([], [], 'b')
def animate(i):
    line.set_data(x, T[i,:])
    return line,


anim = animation.FuncAnimation(fig, animate, frames=t, interval=10, blit=True, repeat=True)

plt.grid()
plt.title('Wave propagation')
plt.xlabel('Distance')
plt.ylabel('Amplitude')

plt.show()
