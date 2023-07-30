"""
Code for the plots of the solutions of the transport equation
The code takes the file with the simulation's parameters
and the solution's files as input from the shell
"""
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

description = 'code for plotting and analyzing the solution'

parser = argparse.ArgumentParser(description=description)
parser.add_argument('parameter', help='path of the file with the input parameters')
parser.add_argument('solutions', help='path of the file with the solutions')
args = parser.parse_args()

# Parameters
par = np.loadtxt(args.parameter)
N   = int(par[0]) + 1
T   = int(par[1]) + 1
v   = par[2]
dt  = par[3]
dx  = par[4]

# load solutions
u = np.loadtxt(args.solutions)
# reshape for plot
sol = np.reshape(u, (T, N))
# x axis
x = np.linspace(0, N*dx, N)

#==========================================================
# Plot
#==========================================================

fig = plt.figure(1)
ax = fig.add_subplot(projection='3d')

ax.set_xlabel('Distance')
ax.set_ylabel('Time')
ax.set_zlabel('Amplitude')


gridx, gridy = np.meshgrid(x, range(T))
ax.plot_surface(gridx, gridy, sol)

plt.figure(2)
freq = np.fft.rfftfreq(N, dx)

for i in range(4):
    plt.subplot(2,2,i+1)
    plt.title(f'Fourier trasform at time {i*(T//4)}')
    sp = np.fft.rfft(sol[i*(T//4),:])
    plt.xlabel("k", fontsize=10)
    plt.ylabel("$|sol(k)|^2$", fontsize=10)
    plt.plot(freq, abs(sp)**2, marker='.', linestyle='-', color='b')
    plt.yscale('log')
    plt.grid()
    
fig = plt.figure(3)

plt.title('Animation of solution', fontsize=15)
plt.xlabel('Distance')
plt.ylabel('Amplitude')
plt.grid()
plt.xlim(np.min(x), np.max(x))
plt.ylim(np.min(sol[0,:]) - 0.3, np.max(sol[0,:]) + 0.3)

dot, = plt.plot([], [], 'b-')

def animate(i):
    
    dot.set_data(x, sol[i,:])
    return dot,

anim = animation.FuncAnimation(fig, animate, frames=np.arange(0, T, 1) ,interval=10, blit=True, repeat=True)

plt.show()
