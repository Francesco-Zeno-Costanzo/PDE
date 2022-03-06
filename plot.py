import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

par = np.loadtxt(r'input_lw.txt')
N = int(par[0]) + 1
T = int(par[1]) + 1
v = par[2]
dt = par[3]
dx = par[4]

u = np.loadtxt(r'tra_lw1.dat')

sol = np.reshape(u, (T, N))

x = np.linspace(0, N*dx, N)

fig = plt.figure(1)
ax = fig.add_subplot(projection='3d')

ax.set_xlabel('Distanza')
ax.set_ylabel('Tempo')
ax.set_zlabel('Ampiezza')


gridx, gridy = np.meshgrid(x, range(T))
ax.plot_surface(gridx, gridy, sol)

plt.figure(2)
freq = np.fft.rfftfreq(N, dx)

for i in range(4):
    plt.subplot(2,2,i+1)
    plt.title(f'spettro soluzione al tempo {i*(T//4)}')
    sp = np.fft.rfft(sol[i*(T//4),:])
    plt.xlabel("k", fontsize=10)
    plt.ylabel("$|sol(k)|^2$", fontsize=10)
    plt.plot(freq, abs(sp)**2, marker='.', linestyle='')
    plt.yscale('log')
    plt.grid()
    
fig = plt.figure(3)

plt.title('Animazione soluzione', fontsize=15)
plt.xlabel('distanza')
plt.ylabel('ampiezza')
plt.grid()
plt.xlim(np.min(x), np.max(x))
plt.ylim(np.min(sol[0,:]) - 0.1, np.max(sol[0,:]) + 0.1)

dot, = plt.plot([], [], 'b-')

def animate(i):
    
    dot.set_data(x, sol[i,:])
    return dot,

anim = animation.FuncAnimation(fig, animate, frames=np.arange(0, T, 1) ,interval=10, blit=True, repeat=True)

plt.show()
