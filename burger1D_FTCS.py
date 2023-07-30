"""
code to solve burger equation using method TFCS
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#=============================================================
# Parameters
#=============================================================

N   = 201
T   = 800
dx  = 2*np.pi/(N-1)
nu  = 0.07
dt  = 0.002
x   = np.linspace(0, 2*np.pi, N)
sol = np.zeros((N, T+1))

# Initial condition
sol[:,0] = np.sin(x)

#=============================================================
# Solutions
#=============================================================

for n in range(T):
    for i in range(1, N - 1):
        # Euler method for second derivative
        ddsol = sol[i+1, n] - 2*sol[i, n] + sol[i-1, n]
        # central difference for first derivative
        dsol = sol[i+1, n] - sol[i-1, n]
        # Eulet fot time evolutions
        sol[i, n+1] = sol[i, n] - sol[i, n]*dt/(2*dx)*dsol + nu*dt/dx**2*ddsol

    # Periodic boundary conditions
    sol[0, n+1] = sol[N-1, n+1]
    sol[N-1, n+1] = sol[1, n+1]

#=============================================================
# Plot
#=============================================================

plt.figure(1)
freq = np.fft.rfftfreq(N, dx/(2*np.pi))
sp_i = np.fft.rfft(sol[:, 0])
sp_f = np.fft.rfft(sol[:,-1])

plt.subplot(211)
plt.title(f'Fourier trasform of solutions at t = 0')
plt.ylabel("$|sol(k)|^2$", fontsize=10)
plt.plot(freq, abs(sp_i)**2, marker='.', linestyle='-', color='b')
plt.ylim(1e-5, 1e5)
plt.yscale('log')
plt.grid()
plt.subplot(212)
plt.title(f'Fourier trasform of solutions at = {T}')
plt.xlabel("k", fontsize=10)
plt.ylabel("$|sol(k)|^2$", fontsize=10)
plt.plot(freq, abs(sp_f)**2, marker='.', linestyle='-', color='b')
plt.yscale('log')
plt.grid()

#=============================================================
# Animations
#=============================================================

fig = plt.figure(2)
plt.xlim(np.min(x), np.max(x))
plt.ylim(np.min(sol), np.max(sol))

line, = plt.plot([], [], 'b')
def animate(i):
    line.set_data(x, sol[:,i])
    return line,


anim = animation.FuncAnimation(fig, animate, frames=np.arange(0, T, 5), interval=1, blit=True, repeat=True)

plt.grid()
plt.title('burger equation')
plt.xlabel('Distanza')
plt.ylabel('ampiezza')

#anim.save('calore.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

plt.show()
