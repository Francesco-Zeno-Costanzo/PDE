"""
Code for solving the heat equation using implicit FTCS
"""
import numpy as np
import matplotlib as mp
import scipy.sparse as sp
import scipy.sparse.linalg
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#===============================================================
# Parameters 
#===============================================================
N     = 100                    # number of points on x axis
x     = np.linspace(0, N, N)   # x axis
tstep = 5000                   # number of points on time axis
D     = 0.5                    # diffusion coefficent
dx    = 0.01                   # step on x
dt    = 1e-4                   # step on time
r     = D*dt/dx**2             # parameter of equations
#r must be lees than 1/2 otherways the integration fail
T = np.zeros((N,tstep))
# Initial condition
T_0 = 500*np.exp(-((50-x)/20)**2)
#T_0 = 500*np.exp(-((5-x)/20)**2) + 500*np.exp(-((95-x)/20)**2)

print(r)

#===============================================================
# Solutions 
#===============================================================

A = sp.diags([-r, 1+2*r, -r], [-1, 0, 1], shape=(N, N)).toarray()

T[:, 0] = T_0
b       = T_0

for i in range(tstep):

    u = np.linalg.solve(A, b)
    b = u
    T[:, i] = b

#===============================================================
# PLOT 
#===============================================================

fig = plt.figure(1)
ax = fig.add_subplot(projection='3d')
gridx, gridy = np.meshgrid(range(tstep), x)
ax.plot_surface(gridx, gridy, T, cmap=mp.cm.coolwarm,vmax=250,linewidth=0,rstride=2, cstride=100)
ax.set_title('Heat diffussion')
ax.set_xlabel('Time')
ax.set_ylabel('Distance')
ax.set_zlabel('Temperature')

#===============================================================
# Animation 
#===============================================================

fig = plt.figure(2)
plt.xlim(np.min(x), np.max(x))
plt.ylim(np.min(T), np.max(T))

line, = plt.plot([], [], 'b')
def animate(i):
    line.set_data(x, T[:,i])
    return line,


anim = animation.FuncAnimation(fig, animate, frames=tstep, interval=10, blit=True, repeat=True)

plt.grid()
plt.title('Heat diffussion')
plt.xlabel('Distance')
plt.ylabel('Temperature')

#anim.save('calore.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

plt.show()
