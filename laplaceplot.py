import numpy as np
import matplotlib.pyplot as plt

A = np.loadtxt('laplace.dat', unpack=True)
N, M = int(A[0]), int(A[1])

b = A[2:M+2]
u = A[M+ 2:]
x = np.linspace(0, 1, N)

B = np.reshape(b, (N, N))
U = np.reshape(u, (N, N))

#Plot
fig = plt.figure(1)
gridx, gridy = np.meshgrid(x, x)
ax = fig.add_subplot(1,2,1, projection='3d')
#plot potenziale generato dalla distribuzione
ax.plot_surface(gridx, gridy, U, color='orange')
ax.set_title('Soluzione equazione di laplace')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('potenziale')

ax1 = fig.add_subplot(1,2,2, projection='3d')

#plot distribuzione di carica
ax1.plot_surface(gridx, gridy, B)
ax1.set_title('condizioni al contorno e sorgente')
ax1.set_xlabel('X')
ax1.set_ylabel('Y')
ax1.set_zlabel('densità')

plt.show()
