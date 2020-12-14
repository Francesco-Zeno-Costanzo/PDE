import numpy as np
import matplotlib as mp
import matplotlib.pyplot as plt
import matplotlib.animation as animation

N=125
t=314
c=1
y=np.linspace(0, N, N)
U=np.zeros((N,t))
U[0:N,0]=np.sin(0.05*y)

for j in np.arange(0,t-2):
    for i in np.arange(0,N-1):
        U[i, j+1] = 2*U[i,j]-U[i,j-1]+ c*(U[i+1, j]-2*U[i, j]+U[i-1, j])


fig = plt.figure()
plt.xlim(np.min(y), np.max(y))
plt.ylim(np.min(U)-0.1, np.max(U)+0.1)

line, = plt.plot([], [], 'b')
def animate(i):
    line.set_data(y, U[:,i])
    return line,


anim = animation.FuncAnimation(fig, animate, frames=t ,interval=20, blit=True, repeat=True)

plt.grid()
plt.title('Propagazione onde')
plt.xlabel('Distanza')
plt.ylabel('Ampiezza')

#anim.save('onda.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

plt.show()