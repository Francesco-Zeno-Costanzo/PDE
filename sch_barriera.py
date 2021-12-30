import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg
import matplotlib.pyplot as plt
import matplotlib.animation as animation


def psi_inc(x):
    '''Funzione d'onda "particella" incidente
    '''

    x0 = -0.1      #centro pacchetto
    a = 0.01       #larghezza pacchetto
    k = 200000.0   #numero donda del pacchetto

    A = 1. #/ np.sqrt( 2 * np.pi * a**2 ) #per avere pacchetto normalizzato
    K1 = np.exp( - ( x - x0 )**2 / ( 2. * a**2 ) )
    K2 = np.exp( 1j * k * x )

    return A * K1 * K2

def Potenziale(x):
    '''Potenziale a "delta"
    '''
    A = 1e6
    s = 0.001
    return A*np.exp(-(x/s)**2)


n = 1001
dx = 0.001
a = -0.5
b = -a
x = np.linspace(a, b, n)


T = 0.0035
dt = 1e-5
t = 0
time_steps = int(T/dt)

k1 = - 1j / (2.)
k2 =  1j


psi = psi_inc(x)

#matrice cinetica
P = (k1 / dx**2) * sp.diags([1, -2, 1], [-1, 0, 1], shape=(n, n))

#matrice potenziale
V = k2 * sp.diags(Potenziale(x), 0, shape=(n, n))

#identità
I =  sp.diags([1], 0, shape=(n, n))

PSI = np.zeros((time_steps, len(psi)))
PSI[0,:] = abs(psi)**2

fig = plt.figure()
plt.title("Propagazione pacchetto gaussiano")

plt.plot(x, Potenziale(x), label='$V(x)$' )

plt.grid()

plt.ylim(-0.1, np.max(PSI[0,:]))


for i in range(time_steps):


    A = (I - dt/2. * (P + V))
    b = (I + dt/2. * (P + V)) * psi

    psi = sp.linalg.spsolve(A,b)

    t += dt

    PSI[i, :] = abs(psi)**2

line, = plt.plot([], [], 'b', label=r"$|\psi(x, t)|^2$")

def animate(i):
    line.set_data(x, PSI[i, :])
    return line,

plt.legend(loc='best')

anim = animation.FuncAnimation(fig, animate, frames=time_steps, interval=10, blit=True, repeat=True)

#anim.save('tunnel barriera.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

plt.show()