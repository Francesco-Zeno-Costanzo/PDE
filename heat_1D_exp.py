"""
Code for solving the heat equation using explicit FTCS
"""
import numpy as np
import matplotlib as mp
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
T[:, 0] = 500*np.exp(-((50-x)/20)**2)
#T[:, 0] = 500*np.exp(-((5-x)/20)**2) + 500*np.exp(-((95-x)/20)**2)

print(r)

#===============================================================
# Solutions 
#===============================================================

for time in range(1,tstep):
    for i in range(1,N-1):
        T[i,time] = T[i,time-1] + r*(T[i-1,time-1]+T[i+1,time-1]-2*T[i,time-1])
#    T[0,time]=T[1,time] # uncomment to have floating edges
#    T[N-1,time]=T[N-2,time]

#===============================================================
# PLOT 
#===============================================================

fig = plt.figure(1)
ax = fig.add_subplot(projection='3d')
gridx, gridy = np.meshgrid(range(tstep), range(N))
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
