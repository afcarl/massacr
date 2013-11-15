import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
import streamplot as sp
plt.rcParams['contour.negative_linestyle'] = 'solid'

#####################
# LOAD MODEL OUTPUT #
#####################

t = np.loadtxt('t1.txt',delimiter='\n')
x0 = np.loadtxt('x1.txt',delimiter='\n')
y0 = np.loadtxt('y1.txt',delimiter='\n')

x=x0
y=y0
bits = len(x)
x = np.append(x0, np.max(x0)+.001)
y = np.append(y0, np.max(y0)+(y0[-1]-y0[-2]))

xg, yg = np.meshgrid(x[:],y[:])

h0 = np.loadtxt('hMat.txt')
u0= np.loadtxt('uMat1.txt')
v0= np.loadtxt('vMat1.txt')
psi0 = np.loadtxt('psiMat.txt')
rho = np.loadtxt('rho1.txt')
viscosity = 1e-3
permeability = np.loadtxt('permeability1.txt')
permeability = permeability


fig=plt.figure()


i=150

print h0.shape

#######################
# MAKE DATA PLOTTABLE #
#######################

h = h0[i*len(y)-i:((i)*len(y)+len(x))-i-1,:]
h = np.append(h, h[-1:,:], axis=0)
h = np.append(h, h[:,-1:], axis=1)

psi = psi0[i*len(y)-i:((i)*len(y)+len(x))-i-1,:]
psi = np.append(psi, psi[-1:,:], axis=0)
psi = np.append(psi, psi[:,-1:], axis=1)

v = v0[i*len(y)-i:((i)*len(y)+len(x))-i-1,:]
v = np.append(v, v[-1:,:], axis=0)
v = np.append(v, v[:,-1:], axis=1)

u = u0[i*len(y)-i:((i)*len(y)+len(x))-i-1,:]
u = np.append(u, u[-1:,:], axis=0)
u = np.append(u, u[:,-1:], axis=1)

permeability = np.append(permeability, permeability[-1:,:], axis=0)
permeability = np.append(permeability, permeability[:,-1:], axis=1)

rho = np.append(rho, rho[-1:,:], axis=0)
rho = np.append(rho, rho[:,-1:], axis=1)


#vmin0=np.min(h[h.shape[0]*1300.0/3000.0:,:])
#vmax0=np.max(h[h.shape[0]*1300.0/3000.0:,:])

# TRYING SOMETHING

h = h[h.shape[0]*1700.0/3000.0:,:]
psi = psi[psi.shape[0]*1700.0/3000.0:,:]
xg = xg[xg.shape[0]*1700.0/3000.0:,:]
yg = yg[yg.shape[0]*1700.0/3000.0:,:]
u = u[u.shape[0]*1700.0/3000.0:,:]
v = v[v.shape[0]*1700.0/3000.0:,:]

####################
# STREAM FUNCTIONS #
####################

ax1=fig.add_subplot(1,1,1, aspect='equal')

p = plt.contour(xg,yg,h,20,cmap=cm.rainbow)
plt.clabel(p,inline=True,fontsize=12,fontweight='bold')

CS = plt.contour(xg, yg, psi, 40, colors='k',linewidths=np.array([1.4]))
#plt.quiver(xg,yg,u,v)

#plt.clabel(CS,CS.levels[CS.levels>0],inline=False,fmt='>',fontsize=12,fontweight='bold')
#plt.clabel(CS,CS.levels[CS.levels<0],inline=False,fmt='<',fontsize=12,fontweight='bold')

#plt.title("STREAMFUNCTIONS",fontsize=8)

plt.xlim(np.min(x), np.max(x))
#plt.ylim(-1300, np.max(y))

    
plt.subplots_adjust(bottom=.2, left=.1, right=.90, top=0.9, hspace=.3)

#cax = fig.add_axes([0.2, 0.1, 0.6, 0.03])
cax = fig.add_axes([0.2, 0.2, 0.6, 0.03])
cbar = plt.colorbar(p, cax=cax,orientation='horizontal',ticks=[])
cbar.set_label(r'TEMPERATURE [K]',fontsize=8)

plt.savefig('n14.png')
print "yeah!"

#####################
# TOP "HEAT" OUTPUT #
#####################

##ax1=fig.add_subplot(1,1,1)
##
##print x.shape
##print h[-1,:].shape
##
##plt.plot(x,h[-1,:])
##
##plt.savefig('heat.png')
##print "yeah!"

