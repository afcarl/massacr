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
feldspar0 = np.loadtxt('feldsparMat.txt')
glass0 = np.loadtxt('glassMat.txt')


fig=plt.figure()


i=9
cell = 4
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

feldspar = feldspar0[i*len(y0)/cell:i*len(y0)/cell+len(y0)/cell,:]
#feldspar = np.append(feldspar, feldspar[-1:,:], axis=0)
#feldspar = np.append(feldspar, feldspar[:,-1:], axis=1)

glass = glass0[i*len(y)/cell-i:((i)*len(y)/cell+len(x)/cell)-i-1,:]
#glassmin = glass0[3*len(y)/cell-3:((3)*len(y)/cell+len(x)/cell)-3-1,:]
#glass = np.append(glass, glass[-1:,:], axis=0)
#glass = np.append(glass, glass[:,-1:], axis=1)

# SELECT RELEVANT PART OF MODEL DOMAIN

#h = h[h.shape[0]*1700.0/3000.0:,:]
#psi = psi[psi.shape[0]*1700.0/3000.0:,:]
#xg = xg[xg.shape[0]*1700.0/3000.0:,:]
#yg = yg[yg.shape[0]*1700.0/3000.0:,:]
#u = u[u.shape[0]*1700.0/3000.0:,:]
#v = v[v.shape[0]*1700.0/3000.0:,:]
feldspar = feldspar[feldspar.shape[0]*1700.0/3000.0:,:]
#glass = glass[glass.shape[0]*1700.0/3000.0:,:]


####################
# STREAM FUNCTIONS #
####################

ax1=fig.add_subplot(2,1,1, aspect='equal')
levels00 = np.linspace(.00015, np.max(psi), 15)
levels0 = np.linspace(np.min(psi), -.00015, 15)
levels = np.append(levels0,levels00,axis=1)
CS = plt.contour(xg, yg, psi, levels, colors='k',linewidths=np.array([2.0]))
#plt.quiver(xg,yg,u,v)

plt.title("STREAMFUNCTIONS",fontsize=8)

plt.xlim(np.min(x), np.max(x))


#############
# ISOTHERMS #
#############

ax1=fig.add_subplot(2,1,2, aspect='equal')
p = plt.contour(xg,yg,h-273.0,np.arange(0,150,10),
                colors='k',linewidths=np.array([2.0]))
plt.clabel(p,inline=True,fontsize=8,fontweight='bold')

plt.title("ISOTHERMS",fontsize=8)

plt.xlim(np.min(x), np.max(x))


########################
# SINGLE CROSS-SECTION #
########################

##print glass[11,:]
##
##ax1=fig.add_subplot(2,2,3)
##p = plt.plot(glass[11,:])
##plt.xlim(0, len(x)/cell-1)
##plt.title("BASALTIC GLASS (horiz. x-section)",fontsize=8)



##############################
# BASALTIC GLASS DISSOLUTION #
##############################

##print glass0.shape
##print glass.shape
##ax1=fig.add_subplot(2,2,4, aspect='equal')
##p = plt.pcolor(glass,cmap=cm.spectral)
##
##plt.title("AMOUNT OF BASALTIC GLASS [mol]",fontsize=8)
##plt.colorbar(p,orientation='horizontal')
##
##plt.xlim(0, len(x)/cell-1)
##plt.ylim(0, len(y)/cell-1)
##
##plt.subplots_adjust(bottom=.2, left=.1, right=.90, top=0.9, hspace=.3)




# COLORBAR STUFF, MAYBE

#cax = fig.add_axes([0.2, 0.1, 0.6, 0.03])
##cax = fig.add_axes([0.2, 0.2, 0.6, 0.03])
#cbar = plt.colorbar(p, cax=cax,orientation='horizontal')
#cbar.set_label(r'TEMPERATURE [K]',fontsize=8)

plt.savefig('d1220.png')
print "flow field plots"

###################
# BENCHMARK PLOTS #
###################

fig=plt.figure()

# TOP HEAT FLUX
ax1=fig.add_subplot(2,1,1)
plt.plot([0,3000],[.27,.27],'r-')
plt.plot(x,-1.6*(h[-3,:]-h[-4,:])/(y[2]-y[1]))
plt.xlim(0,3000)
plt.ylim(0,1)
plt.xlabel('x [m]',fontsize=8)
plt.ylabel('HEAT FLUX [W/m^2]',fontsize=8)
plt.title('HEAT FLUX',fontsize=8)

# TOP FLUID FLUX
ax1=fig.add_subplot(2,1,2)
plt.plot([0,3000],[0.0,0.0],'r-')
plt.plot(x,-v[-1,:])
plt.xlim(0,3000)
#plt.ylim(-10.0e-14,10.0e-14)
#plt.yticks([-10e-12, -5e-12, 0, 5e-12, 10e-12])
plt.xlabel('x [m]',fontsize=8)
plt.ylabel('FLUID FLUX [m/s]',fontsize=8)
plt.title('FLUID FLUX',fontsize=8)

plt.savefig('benchGraphs0.png')
print "benchmark plots"

print sum(sum(np.sqrt(u**2+v**2)))
