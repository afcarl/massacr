import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
import streamplot as sp
plt.rcParams['contour.negative_linestyle'] = 'solid'
plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)

#####################
# LOAD MODEL OUTPUT #
#####################

t = np.loadtxt('t.txt',delimiter='\n')
x0 = np.loadtxt('x.txt',delimiter='\n')
y0 = np.loadtxt('y.txt',delimiter='\n')

x=x0
y=y0
bits = len(x)
x = np.append(x0, np.max(x0)+.001)
y = np.append(y0, np.max(y0)+(y0[-1]-y0[-2]))

xg, yg = np.meshgrid(x[:],y[:])

h0 = np.loadtxt('hMat.txt') # no 1
u0= np.loadtxt('uMat.txt')
v0= np.loadtxt('vMat.txt')
psi0 = np.loadtxt('psiMat.txt') # no 1
feldspar0 = np.loadtxt('feldsparMat.txt') 
glass0 = np.loadtxt('caMat.txt')
perm0 = np.loadtxt('permeability.txt')



    
i=9
cell = 6
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


perm = np.append(perm0, perm0[-1:,:], axis=0)
perm = np.append(perm, perm[:,-1:], axis=1)

v = v0[i*len(y)-i:((i)*len(y)+len(x))-i-1,:]
v = np.append(v, v[-1:,:], axis=0)
v = np.append(v, v[:,-1:], axis=1)


u = u0[i*len(y)-i:((i)*len(y)+len(x))-i-1,:]
u = np.append(u, u[-1:,:], axis=0)
u = np.append(u, u[:,-1:], axis=1)

feldspar = feldspar0[(i*len(y0)/cell):(i*len(y0)/cell+len(y0)/cell)-1,:]
feldspar = np.append(feldspar, feldspar[-1:,:], axis=0)
feldspar = np.append(feldspar, feldspar[:,-1:], axis=1)

glass = glass0[(i*len(y0)/cell):(i*len(y0)/cell+len(y0)/cell)-1,:]
glass = np.append(glass, glass[-1:,:], axis=0)
glass = np.append(glass, glass[:,-1:], axis=1)


####################
# STREAM FUNCTIONS #
####################

fig=plt.figure()

ax1=fig.add_subplot(1,1,1, aspect='equal')
levels00 = np.linspace(.000002, np.max(psi), 15)
levels0 = np.linspace(np.min(psi), -.000002, 15)
levels = np.append(levels0,levels00,axis=1)

# permeability plot
permC = plt.contour(xg, yg, np.log10(perm), [-14.0,-14.1], colors='w',linewidths=np.array([2.0]))
#permC = plt.contourf(xg, yg, np.log10(perm), 10, cmap=cm.summer)

# stream function plot
# levels[::2],
#plt.clabel(CS,  inline=0, fmt='>', fontsize=14)

CS = plt.contour(xg, yg, psi, 10, colors='k',linewidths=np.array([1.0]))


p = plt.contourf(xg,yg,h-272.0, np.arange(0.0,126.0,5.0), cmap=cm.rainbow)
plt.clim(0.0,126.0)
cbar = plt.colorbar(p, orientation='horizontal', ticks=np.arange(0.0,126.0,25.0))
cbar.ax.set_xlabel('FLUID TEMPERATURE [$^{\circ}$C]')


#np.putmask(u, np.abs(u) <= 1.0e-9, 0)
#np.putmask(v, np.abs(v) <= 1.0e-9, 0)
#CS = sp.streamplot(ax1, x, y, u, v, color='k', linewidth=1.0)

#plt.quiver(xg,yg,u,v)
plt.yticks([0.0, -500.0, -1000.0], [0, -500, -1000.0])
plt.xticks([0.0, 1500.0, 3000.0], [0, 1500, 3000])

#plt.title("STREAMFUNCTIONS",fontsize=8)

plt.xlim(np.min(x), np.max(x))

plt.savefig('j14.png')


####################
# GEOCHEM CONTOURS #
####################

fig=plt.figure()

ax1=fig.add_subplot(1,1,1, aspect='equal')

# glass plot
xCell = x0
yCell = y0
xCell = xCell[::cell]
yCell= yCell[::cell]

print xCell
xCell = np.append(xCell, np.max(xCell)+.001)
print xCell
yCell = np.append(yCell, np.max(yCell)+.001)
print yCell


pGlass = plt.contourf(xCell, yCell[:-1],glass, 20, cmap=cm.rainbow)

#pGlass = plt.contourf(xg, yg, v, 20, cmap=cm.rainbow)

cbar= plt.colorbar(pGlass, orientation='horizontal')
cbar.ax.set_xlabel('AMOUNT OF BASALTIC GLASS [mol]')


#plt.savefig('expCapNextNext0'+str(i)+'.png')

plt.savefig('j15.png')


print "flow field plots"


###################
# BENCHMARK PLOTS #
###################



fig=plt.figure()


# TOP FLUID FLUX
ax1=fig.add_subplot(2,1,1)
plt.plot([0,np.max(x)],[0.0,0.0],'k--')
p1 = plt.plot(x,-v[-1,:],'gold',linewidth=3, label='SEPARATED OUTCROPS')
plt.xlim(0,np.max(x))
#plt.ylim(-10.0e-12,10.0e-12)
#plt.yticks([-10e-12, -5e-12, 0, 5e-12, 10e-12],[-10e-12, -5e-12, 0, 5e-12, 10e-12])
plt.ylabel('FLUID FLUX [NORMALIZED]',fontsize=10)
plt.xlabel('x DISTANCE [m]',fontsize=10)
handles, labels = ax1.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.legend(handles, labels,loc='best',prop={'size':8}, ncol=1)
plt.savefig('expoFlux0'+str(i)+'.png')

print "benchmark plots"

print sum(sum(np.sqrt(u**2+v**2)))
