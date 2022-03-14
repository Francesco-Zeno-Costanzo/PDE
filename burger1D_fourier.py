import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import matplotlib.animation as animation


L = 1
T = 0.31

dx = 0.001
dt = 0.001
nu = 0.005

Nx = int(L/dx)
Nt = int(T/dt)

t = np.linspace(0, T, Nt)
x = np.linspace(0, L, Nx)

#vettori d'onda
k = 2*np.pi*np.fft.rfftfreq(Nx, dx)
#condizione inziale
u0 = np.sin(2*np.pi*x)

def eq(u, t, k, nu):
    '''
    equazione da risolvere in trasformata spaziale
    '''
    #passo in trasformata così le
    #derivare sono moltiplicazioni
    u_hat = np.fft.rfft(u)
    du_hat = 1j*k*u_hat
    ddu_hat = -k**2*u_hat

    #antitrasformo per avere la derivta
    du = np.fft.irfft(du_hat)
    ddu = np.fft.irfft(ddu_hat)

    #adesso la pde è diventata una ode
    u_t = -u*du + nu*ddu

    return u_t.real


sol = odeint(eq, u0, t, args=(k, nu,)).T

##plot
plt.figure(1)
freq = np.fft.rfftfreq(Nx, dx)
sp_i = np.fft.rfft(sol[:, 0])
sp_f = np.fft.rfft(sol[:,-1])

plt.subplot(211)
plt.title(f'spettro soluzione a t = 0')
plt.ylabel("$|sol(k)|^2$", fontsize=10)
plt.plot(freq, abs(sp_i)**2, marker='.', linestyle='')

plt.yscale('log')
plt.grid()
plt.subplot(212)
plt.title(f'spettro soluzione a t = {T}')
plt.xlabel("k", fontsize=10)
plt.ylabel("$|sol(k)|^2$", fontsize=10)
plt.plot(freq, abs(sp_f)**2, marker='.', linestyle='')
plt.yscale('log')
plt.grid()

fig = plt.figure(2)
plt.xlim(np.min(x), np.max(x))
plt.ylim(np.min(sol), np.max(sol))

line, = plt.plot([], [], 'b')
def animate(i):
    line.set_data(x, sol[:,i])
    return line,


anim = animation.FuncAnimation(fig, animate, frames=Nt, interval=5, blit=True, repeat=True)

plt.grid()
plt.title('burger equation')
plt.xlabel('Distanza')
plt.ylabel('ampiezza')

#anim.save('buger.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

plt.show()


