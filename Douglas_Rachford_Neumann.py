# -*- coding: utf-8 -*-
#Author: Diogo Rocha
#Year: 2016


# -*- coding: utf-8 -*-


from scipy.integrate import simps
from scipy.sparse.linalg import inv
from scipy.sparse.linalg import spsolve
import numpy as np
import scipy as sp
import pylab as pl
from scipy import ndimage
import scipy.sparse as sps
from mpl_toolkits.mplot3d import Axes3D
#from mayavi.mlab import *



#FUNCAO DE ONDA
def phi_inicial(x,y,z):

    #boundary conditions
    #norma = np.real(1.0/simps(simps(gauss,axis=0)))
    #return norma*gauss
    wavefunction = np.sin(np.pi*x)*np.sin(np.pi*y)*np.sin(np.pi*z)

    return wavefunction

def phi_theoretical(x,y,z,t):
    theowave = np.exp(-3*t*miu*np.pi**2)*np.sin(np.pi*x)*np.sin(np.pi*y)*np.sin(np.pi*z)
    return theowave

#Esfera dentro de um cubo
def esfera(x,y,z):
    fn=np.zeros((size_box,size_box,size_box))
    fn[np.where((x-x0)**2 +(y-y0)**2 + (z-z0)**2 < R**2)]= const_esf
    return fn

#tamanho da caixa
size_box = 24
L=1
R=0.5
x0=0.5
y0=0.5
z0=0.5
const_esf = 10

#Constante de difusao
#miu = 0.001
miu=0.002
#MESHGRID
x=np.linspace(0,L,size_box)
y=np.linspace(0,L,size_box)
z=np.linspace(0,L,size_box)
X,Y,Z = np.meshgrid(x,y,z)

dx = X[0,1,0]-X[0,0,0]
dy = Y[1,0,0]-Y[0,0,0]
dz = Z[0,0,1]-Z[0,0,0]



#DEFINIR O TEMPO
dt = 5
tfinal = 1.0e3 * dt
time = np.arange(0.0,tfinal,dt)

#Laplaciano
f=np.zeros((len(x),len(y),len(z)))
def Ax(funcao):
    f == funcao
    novafuncao = np.roll(f,-1,2)+ np.roll(f,1,2) - 2*f
    #       aresta da esquerda
#
#
#       o-o-o     1  -2  1
#       o
    for k in range(0,size_box):
        for m in range(0,size_box):

                #o-o-o canto esquerdo


            novafuncao[k,m,0]=f[k,m,0]-2*f[k,m,1]+f[k,m,2]

# aresta da direita

    for k in range(0,size_box):
        for m in range(0,size_box):
            novafuncao[k,m,size_box-1]= f[k,m,size_box-1]-2*f[k,m,size_box-2]+f[k,m,size_box-3]
# aresta de baixo

    for k in range(0,size_box):
        for m in range(1,size_box-1):
            novafuncao[k,size_box-1,m]= f[k,size_box-1,m-1]-2*f[k,size_box-1,m]+f[k,size_box-1,m+1]

#aresta de cima

    for k in range(0,size_box):
        for m in range(1,size_box-1):
            novafuncao[k,0,m] = f[k,0,m-1]-2*f[k,0,m]+f[k,0,m+1]











    return novafuncao

def Ay(funcao):
    f == funcao
    novafuncao = np.roll(f,-1,1)+ np.roll(f,1,1) - 2*f
    #       aresta da esquerda
#
#
#       o-o-o     1  -2  1
#       o
    for k in range(0,size_box):
        for m in range(1,size_box-2):

                #o-o-o canto esquerdo


            novafuncao[k,m,0]=f[k,m,0]-2*f[k,m+1,0]+f[k,m+2,0]

# aresta da direita

    for k in range(0,size_box):
        for m in range(1,size_box-2):
            novafuncao[k,m,size_box-1]= f[k,m,size_box-1]-2*f[k,m+1,size_box-1]+f[k,m+2,size_box-1]
# aresta de baixo

    for k in range(0,size_box):
        for m in range(0,size_box):
            novafuncao[k,size_box-1,m]= f[k,size_box-1,m]-2*f[k,size_box-2,m]+f[k,size_box-3,m]

#aresta de cima

    for k in range(0,size_box):
        for m in range(0,size_box):
            novafuncao[k,0,m] = f[k,0,m]-2*f[k,1,m]+f[k,2,m]

    return novafuncao


def Az(funcao):
    f == funcao
    novafuncao = np.roll(f,-1,0)+ np.roll(f,1,0) - 2*f
    #       aresta da esquerda
#
#
#       o-o-o     1  -2  1
#       o
    for k in range(1,size_box-1):
        for m in range(0,size_box):

                #o-o-o canto esquerdo


            novafuncao[k,m,0]=f[k-1,m,0]-2*f[k,m,0]+f[k+1,m,0]

# aresta da direita

    for k in range(1,size_box-1):
        for m in range(1,size_box-1):
            novafuncao[k,m,size_box-1]= f[k-1,m,size_box-1]-2*f[k,m,size_box-1]+f[k+1,m,size_box-1]
# aresta de baixo

    for k in range(1,size_box-1):
        for m in range(0,size_box):
            novafuncao[k,size_box-1,m]= f[k-1,size_box-1,m]-2*f[k,size_box-1,m]+f[k+1,size_box-1,m]

#aresta de cima

    for k in range(1,size_box-1):
        for m in range(0,size_box):
            novafuncao[k,0,m] = f[k-1,0,m]-2*f[k,0,m]+f[k+1,0,m]

    return novafuncao









#sistema ax=b  2 membro
F = np.zeros((len(x),len(y),len(z)))
fz = np.zeros((len(x),len(y),len(z)))
I = np.zeros((len(x),len(y),len(z)))

#identity 3 dim case
#for q in range(len(size_box)):
	#i[q,q,q]=1.0

q = esfera(X,Y,Z)

onda1= np.zeros((len(x),len(y),len(z)))
onda2= np.zeros((len(x),len(y),len(z)))
onda3= np.zeros((len(x),len(y),len(z)))
newphi= np.zeros((len(x),len(y),len(z)))
phitheoretical = np.zeros((len(x),len(y),len(z)))
fh = np.zeros((len(x),len(y),len(z)))


Identidade = np.identity(size_box**2)

#lxly = -4.0*np.eye(size_box**2) + np.eye(size_box**2,k= 1) + np.eye(size_box**2,k= -1)
lx = -2.0*np.eye(size_box**2) + np.eye(size_box**2,k= 1) + np.eye(size_box**2,k= -1)
lx[0,0]=-1
lx[0,1]=-1.0
lx[0,2]=1
lx[size_box**2-1,size_box**2-1]=1
lx[size_box**2-1,size_box**2-2]=-2


ly = -2.0*np.eye(size_box**2) + np.eye(size_box**2,k= 1) + np.eye(size_box**2,k= -1)
ly[0,0]=1.0
ly[0,1]=-2.0
ly[0,2]=1
ly[size_box**2-1,size_box**2-1]=1
ly[size_box**2-1,size_box**2-2]=-2



pn = 5

F = phi_inicial(X,Y,Z)

#boundary conditions

a = 0.2
F[0] = F[1]-a*dz
F[:,0,:]=F[:,1,:]-a*dy
F[:,:,0]= F[:,:,1]-a*dx

F[size_box-1]=F[size_box-2]-a*dz
F[:,size_box-1,:]=F[:,size_box-2,:]-a*dy
F[:,:,size_box-1]=F[:,:,size_box-2] - a*dx


fstream=open('/Users/diogorocha/Documents/MCE/Projeto/centrevalues/centreneumann.txt','a')
centre=np.zeros(100)

iteo=0
for t in range(1,102):
    phitheoretical = phi_theoretical(X,Y,Z,t)
    phitheoretical[0] = a
    phitheoretical[:,0,:]=a
    phitheoretical[:,:,0]=a
    phitheoretical[size_box-1]=a
    phitheoretical[:,size_box-1,:]=a
    phitheoretical[:,:,size_box-1]=a






    #Douglas-Rachford Algorithm for Three-Dimensions Implementation
    #define pn as a constant


    #b1 is the second term of first step
    # axis convention
    # axis 1 = x
    # axis 2 = y
    # axis 0 = z

    b1 = (-(miu/dx**2)*Ax(F)
    -2*(miu/dy**2)*Ay(F)
    -2*(miu/dz**2)*Az(F)
    + 2*q
    + pn*F)



#  ax=b system discover onda1
    for i in range(size_box):
        c = np.reshape(b1[i],(size_box**2,1))
        A = lx + pn*Identidade
        G = np.linalg.solve(A,c)
        onda1[i] = np.reshape(G, ((size_box,size_box)))

# b2 is the second term of second step
    b2 = (-(miu/dx**2)*Ax(F)
    -(miu/dy**2)*Ay(F)
    -2*(miu/dz**2)*Az(F)

    + pn*F
    - (miu/dx**2)*Ax(onda1)
          + 2*q )


#ax = b systema discover onda2
    for i in range(size_box):
        s = np.reshape(b2[i],(size_box**2,1))
        A = ly + pn*Identidade
        d = np.linalg.solve(A,s)
        onda2[i] = np.reshape(d, ((size_box,size_box)))

# b3 is the second term of third step

    b3= (-(miu/dx**2)*Ax(F)
    -(miu/dy**2)*Ay(F)
    -(miu/dz**2)*Az(F)

    + pn*F
    - (miu/dx**2)*Ax(onda1)
    - (miu/dy**2)*Ay(onda2)
    + 2*q)




# Rotation of the tridimensional matrix
    for g in range(size_box):
        for h in range(size_box):
            fz[h,g] = b3[g,h]

# ax=b system discover onda3
    for i in range(size_box):
        w = np.reshape(fz[i],(size_box**2,1))
        A = lx + pn*Identidade
        d = np.linalg.solve(A,w)
        onda3[i] = np.reshape(d, ((size_box,size_box)))

    # Rotation of the tridimensional matrix
    for g in range(size_box):
        for h in range(size_box):
            fh[g,h] = onda3[h,g]

    F == fh
    pn = pn +5




    #points3d(X, Y, Z, onda3, colormap="copper", scale_factor=.07)




    fig = pl.figure()
    if iteo<100:
        centre[iteo]= fh[size_box/2,size_box/2,size_box/2]
        print('Representa' , iteo,'figuras')
        #ax = fig.add_subplot(111,projection='3d')
    #       plot_trisurf
        pl.clf()
        #ax.plot_surface(Y[size_box/2],Z[size_box/2],onda3[size_box/2],color='blue')
        pl.pcolormesh(X[size_box/2],Z[size_box/2],fh[size_box/2])
    ##
        pl.xlabel('X')
        pl.ylabel('Y')
        pl.colorbar()
        fname='/Users/diogorocha/Desktop/relatório/Discretoneu/DOUGLASRACHFORDneumann%04d.png'%(t+1)
        pl.savefig(fname)
        iteo += 1
        
np.savetxt(fstream,centre)
fstream.close() 









#Laplaciano tridimensional
#lz = np.zeros((3,size_box,size_box))
#a[0,1,1]=1
#a[1,1,0]=1
#a[1,0,1]=1
#a[1,1,2]=1
#a[1,2,1]=1
#a[1,1,1]=-6
#a[2,1,1]=1

#lapl= sps.kron (primeira,sps.identity(size_box))
#primeira = sps.diags ([0,1,0],[-1,0,1],shape=(size_box,size_box),format='csc')
#segunda = sps.diags([2,-6,2],[-1,0,1],shape=(size_box,size_box),format='csc')
#terceira = sps.diags ([0,1,0],[-1,0,1],shape=(size_box,size_box),format='csc')
#lz[0] = np.eye(size_box)
#lz[1] = -2.0*np.eye(size_box)
#lz[2] = np.eye(size_box)












