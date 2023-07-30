"""
code for solving the wave equation in one spatial dimension using the FTCS method
"""
import numpy as np
import matplotlib as mp
import matplotlib.pyplot as plt
import matplotlib.animation as animation


N = 125  # point on spatial  axis
t = 315  # point on temporal axis

#c must be less equal than 1
#c = v*dt**2/dx**2
c = 1
y = np.linspace(0, 1, N)
U = np.zeros((N,t))
U[0:N,0] = np.sin(2*np.pi*y) # initial condition

for j in np.arange(0,t-1):
    for i in np.arange(0,N-1):
        U[i, j+1] = 2*U[i,j]-U[i,j-1]+ c*(U[i+1, j] - 2*U[i, j] + U[i-1, j])

#===============================================================================
# The simulation returns a wave whose amplitude is greater than initial
# condition. This is probably due to problems with the method used
#===============================================================================

fig = plt.figure(1)
ax = fig.add_subplot(projection='3d')

ax.set_title('Wave propagation')
ax.set_xlabel('Time')
ax.set_ylabel('Distance')
ax.set_zlabel('Amplitudde')


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
plt.title('Wave propagation')
plt.xlabel('Distance')
plt.ylabel('Amplitudde')

#anim.save('onda.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

plt.show()
