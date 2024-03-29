import numpy as np
import matplotlib.pyplot as plt

A = np.loadtxt('lap_c.txt', unpack=True)
N = int(A[0]) + 1
u = A[1:]
x = np.linspace(0, 1, N)
U = np.reshape(u, (N, N))

# Compute of elettrical field
Ex = np.zeros(U.shape)
Ey = np.zeros(U.shape)
dx = 1/N

for i in range(N):
    for j in range(N):
        if i >= 1 and  i <N-1:
            Ey[i, j] = -(U[i+1, j] - U[i-1, j])/(2*dx)
        elif i >= 1:
            Ey[i, j] = -(U[i, j] - U[i-1, j])/dx
        else:
            Ey[i, j] = -(U[i+1, j] - U[i, j])/dx
            
        if j >= 1 and j < N-1:
            Ex[i, j] = -(U[i, j+1] - U[i, j-1])/(2*dx)
        elif j>=1:
            Ex[i, j] = -(U[i, j] - U[i, j-1])/dx
        else:
            Ex[i, j] = -(U[i, j+1] - U[i, j])/dx
        

#Plot
fig = plt.figure(1)
gridx, gridy = np.meshgrid(x, x)
ax = fig.add_subplot(1,1,1, projection='3d')
#plot potenziale generato dalla distribuzione
ax.plot_surface(gridx, gridy, U, color='orange')
ax.set_title('Solution of laplace equations')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('potential')

plt.figure(2)
plt.contour(gridx, gridy, U, cmap='plasma')
plt.streamplot(gridx, gridy, Ex, Ey, density=2)
plt.xlabel('x')
plt.ylabel('y')
plt.title('contur of potential and elettrical field')
plt.show()
