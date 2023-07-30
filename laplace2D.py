import numpy as np
import matplotlib.pyplot as plt

# Number of grid points, spacing on x and y is the same
N = 50
x = np.linspace(0, 1, N)

# Matrix created according to the lexicographic ordering
M = N**2
A = np.zeros((M, M))

for i in range(M):
    A[i, i] = -4

for i in range(M-1) :
    A[i, i+1] = 1
    A[i+1, i] = 1

for i in range(1, N):
    A[i*N-1, i*N] = 0
    A[i*N, i*N-1] = 0

for i in range(M-N):
    A[i, i+N] = 1
    A[i+N, i] = 1

#boundary conditions
#and source function
B = np.zeros((N,N))

#souce
y = np.zeros(N)
y[N//4]= 1
y[3*N//4] = -1
B[:,:] = y#np.sin(2*np.pi*x)
B[:,:] = B[:,:]*np.transpose(B[:,:])

#boundary
B[0,:] = 0#np.sin(2*np.pi*x)
B[:,0] = 0#np.sin(2*np.pi*x)
B[:,-1] = 0#np.sin(2*np.pi*x)
B[-1,:] = 0#np.sin(2*np.pi*x)

#transform matrix into array
b = np.reshape(B, M)
#solve
u = np.linalg.solve(A, b)
#back to matrix
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
