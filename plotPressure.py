import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
plt.rcParams['contour.negative_linestyle'] = 'solid'

t = np.loadtxt('t1.txt',delimiter='\n')
x0 = np.loadtxt('x1.txt',delimiter='\n')
#x0=x0/np.max(x0)
y0 = np.loadtxt('y1.txt',delimiter='\n')
#y0=y0/np.max(np.abs(y0))

x=x0
y=y0
bits = len(x)
x = np.append(x0, np.max(x0)+.001)
y = np.append(y0, np.max(y0)+(y0[-1]-y0[-2]))


xg, yg = np.meshgrid(x[:],y[:])


h = np.loadtxt('h1.txt')
u= np.loadtxt('uMat1.txt')
v= np.loadtxt('vMat1.txt')
psi = np.loadtxt('psiMat1.txt')
#pr = np.loadtxt('pMat.txt')
rho = np.loadtxt('rho1.txt')
viscosity = 1e-3
permeability = np.loadtxt('permeability1.txt')
permeability = permeability

i=3
#wut = u0[i*len(y):((i)*len(y)+len(x)),:]

fig=plt.figure()

count = 1
i=99
#, aspect='equal'
ax1=fig.add_subplot(1,1,1, aspect=1.0)


##h = u0[i*len(y)-i:((i)*len(y)+len(x))-i-1,:]
##h = np.append(h, h[-1:,:], axis=0)
##h = np.append(h, h[:,-1:], axis=1)
##
##psi0 = psi[i*len(y)-i:((i)*len(y)+len(x))-i-1,:]
##psi0 = np.append(psi0, psi0[-1:,:], axis=0)
##psi0 = np.append(psi0, psi0[:,-1:], axis=1)
##
##v = v[i*len(y)-i:((i)*len(y)+len(x))-i-1,:]
##v = np.append(v, v[-1:,:], axis=0)
##v = np.append(v, v[:,-1:], axis=1)
##
##u = u[i*len(y)-i:((i)*len(y)+len(x))-i-1,:]
##u = np.append(u, u[-1:,:], axis=0)
##u = np.append(u, u[:,-1:], axis=1)

h = np.append(h, h[-1:,:], axis=0)
h = np.append(h, h[:,-1:], axis=1)

psi = np.append(psi, psi[-1:,:], axis=0)
psi = np.append(psi, psi[:,-1:], axis=1)

v = np.append(v, v[-1:,:], axis=0)
v = np.append(v, v[:,-1:], axis=1)

u = np.append(u, u[-1:,:], axis=0)
u = np.append(u, u[:,-1:], axis=1)

#porosity = np.append(porosity, porosity[-1:,:], axis=0)
#porosity = np.append(porosity, porosity[:,-1:], axis=1)

permeability = np.append(permeability, permeability[-1:,:], axis=0)
permeability = np.append(permeability, permeability[:,-1:], axis=1)

rho = np.append(rho, rho[-1:,:], axis=0)
rho = np.append(rho, rho[:,-1:], axis=1)

p = plt.pcolor(xg,yg,np.log10(permeability),cmap=cm.summer)
#p = plt.pcolor(xg,yg,np.gradient(h)[1],cmap=cm.summer)
print np.gradient(h)[1].shape
#q = plt.contour(xg,yg,np.log10(permeability),cmap=cm.Spectral,linewidths=np.array([3.0]))
#plt.rcParams['contour.negative_linestyle'] = 'dashed'
#p = plt.pcolor(xg,yg,np.log10(permeability+1e-20),cmap=cm.Spectral)
#plt.clim(-17,-12)
print xg.shape
CS = plt.contour(xg, yg, psi, 50, colors='k',linewidths=np.array([1.6]))
CS = plt.contour(xg, yg, h, 35, colors='#660066',linewidths=np.array([2.0]))
plt.clabel(CS, fontsize=9, inline=1,fmt='%3.0f')
#CS = plt.contour(xg, yg, rho, 20, colors='m',linewidths=np.array([1.0]))

#plt.axis('scaled')
#plt.scatter(xv,yv)
#sp.streamplot(ax1,xg, yg, (u),(v))


N = np.ones(u.shape)
for i in range(len(xg)):
    for j in range(len(yg)):
        if v[i,j] != 0.0 or u[i,j]!=0.0:
            N[i,j] = np.sqrt(u[i,j]**2+v[i,j]**2)
#N = np.sqrt(u**2+v**2)

#u=u/N
#v=v/N
widths = np.linspace(0, 2, x.size)
#plt.quiver(xg[::1], yg[::1], 5.0*u[::1],5.0*v[::1],linewidths=[.3],pivot='mid')

plt.title("",fontsize=8)

count = count + 1
plt.xlim(np.min(x), np.max(x))
plt.ylim(-1300, np.max(y))
#plt.ylim(np.min(y), np.max(y))


    
plt.subplots_adjust(bottom=.2, left=.1, right=.90, top=0.9, hspace=.3)

cax = fig.add_axes([0.2, 0.1, 0.6, 0.03])
cbar = plt.colorbar(p, cax=cax,orientation='horizontal')

uniques = np.unique(np.log10(permeability))
print len(uniques)
cbar = plt.colorbar(p, cax=cax,orientation='horizontal')
cbar.set_label(r'log of permeability',fontsize=8)

plt.savefig('sept16.png')
print "yeah!"



