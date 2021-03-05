import time
import glob
import numpy as np
import imageio as io
from PIL import Image
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

start_time=time.time()

N=100
t=600
c=0.5
y=np.linspace(0, N, N)
U=np.zeros((N,N,t))
U[:,:,0]=np.exp(-((y-N//2)**2)/100)
U[:,:,0]=-U[:,:,0]*np.transpose(U[:,:,0])
#per poter vedere la configurazione iniziale
'''
fig = plt.figure(1)
ax = fig.gca(projection='3d')
gridx , gridy = np.meshgrid(range(N), range(N))
ax.plot_surface(gridx,gridy,U[:,:,0])
plt.show()
'''
##
for j in np.arange(0,t-1):
    for i in np.arange(0,N-1):
        for k in np.arange(0,N-1):
            U[i,k,j+1] = 2*U[i,k,j]-U[i,k,j-1]+ c*(U[i+1,k,j]-2*U[i,k,j]+U[i-1,k,j]+U[i,k+1,j]-2*U[i,k,j]+U[i,k-1,j])

a=(time.time() - start_time)
print("--- %s secondi per il calcolo ---" %a)
#si crea un plot per ogni instante di tempo e tutti vengono salvati in una cartella
gridx , gridy = np.meshgrid(range(N), range(N))

ZM=np.max(U)
zm=np.min(U)
for h in range(t):
    fig = plt.figure(h)
    ax = fig.gca(projection='3d')
    ax.set_zlim(zm, ZM)
    ax.set_xlabel('Distanza')
    ax.set_ylabel('Distanza')
    ax.set_zlabel('Ampiezza')
    ax.plot_surface(gridx,gridy,U[:,:,h], cmap=cm.coolwarm)
    ax.set_title('Tamburo')
    plt.savefig(r'C:\Users\franc\Desktop\codici python\gif/%d'%(h))
    plt.close(fig)

b=(time.time() - start_time-a)
print("--- %s secondi per i plot ---" %b)

#si prendono dalla cartella i plot e si uniscono creando l'animazione.
#Essa viene salvata nella stessa cartella dove è il codice
frames=[]
imgs=sorted(glob.glob(r'C:\Users\franc\Desktop\codici python\gif/*.png'))
imgs.sort(key=len)
for i in imgs:
    new_frame=Image.open(i)
    frames.append(new_frame)

frames[0].save('animazione1.gif',format='GIF',append_images=frames[:],save_all=True,duration=50,loop=0)

l=(time.time() - start_time-b-a)
print("--- %s secondi per creare l'animazione ---" %l)