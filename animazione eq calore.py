import numpy as np
import matplotlib as mp
import matplotlib.pyplot as plt
import matplotlib.animation as animation
N=100
x=np.linspace(0, N, N)
tstep=10000
T=np.zeros((N,tstep))
#Temperatura iniziale
T[0:N,0]=500*np.exp(-((50-x)/20)**2)

ts=1
tau=1*ts


for time in range(1,tstep):
    for i in range(0,N-1):
        T[i,time]=T[i,time-1]+0.5*tau/ts*(T[i-1,time-1]+T[i+1,time-1]-2*T[i,time-1])
    T[0,time]=T[1,time]
    T[N-1,time]=T[N-2,time]


fig = plt.figure()

redDot, = plt.plot(T[0:,0]+0.5*tau/ts*(T[0:,0]+T[0:,0]-2*T[0:,0]))

def animate(i):
    redDot.set_data(x, T[:,i-1]+0.5*tau/ts*(T[:,i-1]+T[:,i-1]-2*T[:,i-1]))
    return redDot,


anim = animation.FuncAnimation(fig, animate, frames=range(0, tstep),interval=10, blit=True, repeat=True)

plt.grid()
plt.title('Diffusione del calore')
plt.xlabel('Distanza')
plt.ylabel('Temperatura')

#anim.save('basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

plt.show()
