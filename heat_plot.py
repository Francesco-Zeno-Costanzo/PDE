"""
Code for the plots of the solutions of the heat equation
The code takes the file with the simulation's parameters
and the solution's files as input from the shell
"""
import argparse
import numpy as np
import matplotlib as mp
import matplotlib.pyplot as plt
import matplotlib.animation as animation

description = 'code for plotting and analyzing the solution'

parser = argparse.ArgumentParser(description=description)
parser.add_argument('parameter', help='path of the file with the input parameters')
parser.add_argument('solutions', help='path of the file with the solutions')
args = parser.parse_args()

# Parameters
par = np.loadtxt(args.parameter)
N   = int(par[0])
t   = int(par[1]) + 1
D   = par[2]
dt  = par[3]
dx  = par[4]

# load solutions
u = np.loadtxt(args.solutions)
# reshape for plot
T = np.reshape(u, (t, N))
# x axis
x = np.linspace(0, N*dx, N)

#==========================================================
# Plot
#==========================================================

fig = plt.figure(1)
ax = fig.add_subplot(projection='3d')
gridx, gridy = np.meshgrid(x, range(t))
ax.plot_surface(gridx, gridy, T, cmap=mp.cm.coolwarm,vmax=250,linewidth=0,rstride=2, cstride=2)
ax.set_title('Heat diffussion')
ax.set_ylabel('Time')
ax.set_xlabel('Distance')
ax.set_zlabel('Temperature')

#===============================================================
# Animation 
#===============================================================

fig = plt.figure(2)
plt.xlim(np.min(x), np.max(x))
plt.ylim(np.min(T), np.max(T))

line, = plt.plot([], [], 'b')
def animate(i):
    line.set_data(x, T[i,:])
    return line,


anim = animation.FuncAnimation(fig, animate, frames=t, interval=10, blit=True, repeat=True)

plt.grid()
plt.title('Heat diffussion')
plt.xlabel('Distance')
plt.ylabel('Temperature')

plt.show()
