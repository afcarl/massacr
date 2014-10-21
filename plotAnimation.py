import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
import streamplot as sp
plt.rcParams['contour.negative_linestyle'] = 'solid'
plt.rc('font', family='Arial')

plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)
plt.rcParams['xtick.major.size'] = 0
plt.rcParams['ytick.major.size'] = 0
#plt.rcParams['xtick.direction'] = 'out'
#plt.rcParams['ytick.direction'] = 'out'
plt.rcParams['xtick.major.pad'] = 3
plt.rcParams['ytick.major.pad'] = 3
plt.rcParams['axes.linewidth'] = 1.0
plt.rcParams['axes.color_cycle'] = "#CE1836, #F85931, #EDB92E, #A3A948, #009989"

print "doing something..."

#####################
# LOAD MODEL OUTPUT #
#####################

cell = 10
#steps = 400
steps =4

path = "output/big/"
path = ""


#htest = np.loadtxt(path + 'steady_h.txt')

##fig=plt.figure()
##ax1=fig.add_subplot(1,1,1, aspect='equal')
##ht = plt.contourf(htest)
##cbar = plt.colorbar()
##plt.savefig(path + 'htest.png')
##print "yeah!"


# load output
#t = np.loadtxt(path + 't.txt',delimiter='\n')
x0 = np.loadtxt(path + 'x.txt',delimiter='\n')
y0 = np.loadtxt(path + 'y.txt',delimiter='\n')

# format plotting geometry
x=x0
y=y0
bits = len(x)
#x = np.append(x0, np.max(x0)+.001)
#y = np.append(y0, np.max(y0)+(y0[-1]-y0[-2]))

xCell = x0
yCell = y0
xCell = xCell[::cell]
yCell= yCell[::cell]
xCell = np.append(xCell, np.max(xCell)+xCell[0])
yCell = np.append(yCell, np.max(yCell)-yCell[-1])

print xCell
print yCell

xg, yg = np.meshgrid(x[:],y[:])


# load output
#h0 = np.loadtxt(path + 'hMat.txt')
#u0= np.loadtxt(path + 'uMat.txt')
#v0= np.loadtxt(path + 'vMat.txt')
#psi0 = np.loadtxt(path + 'psiMat.txt')
feldspar0 = np.loadtxt(path + 'pri_feldspar.txt') 
glass0 = np.loadtxt(path + 'pri_glass.txt')
perm0 = np.loadtxt(path + 'permeability.txt')
geo0 = np.loadtxt(path + 'sol_c.txt')

# format output
geo00 = np.zeros(steps)
print "go"
for i in range(steps):
    geo1 = geo0[(i*len(y0)/cell):(i*len(y0)/cell+len(y0)/cell)-1,:]
    geo1 = np.append(geo1, geo1[-1:,:], axis=0)
    geo1 = np.append(geo1, geo1[:,-1:], axis=1)
    geo00[i] = np.max(geo1)


# load output
steady_h = np.loadtxt(path + 'steady_h.txt')
steady_u = np.loadtxt(path + 'steady_u.txt')
steady_v = np.loadtxt(path + 'steady_v.txt')
steady_psi = np.loadtxt(path + 'steady_psi.txt')
steady_permeability = np.loadtxt(path + 'steady_permeability.txt')

##############
# FIRST PLOT #
##############

fig=plt.figure()
ax1=fig.add_subplot(1,1,1, aspect='equal')

# plot permeability contours
permC = plt.contour(xg, yg, np.log10(steady_permeability), [-14.0,-14.1], colors='w',linewidths=np.array([2.0]))
#plt.clabel(CS,  inline=0, fmt='>', fontsize=14)

# plot streamfunctions
print xg.shape
print yg.shape
print steady_psi.shape
CS = plt.contour(xg, yg, steady_psi, 10, colors='k',linewidths=np.array([1.0]))

# plot temperature
print xg.shape
print yg.shape
print steady_h.shape
p = plt.contourf(xg,yg,steady_h-272.0, np.arange(0.0,np.max(steady_h-272.0),5.0), cmap=cm.rainbow)
plt.clim(0.0,np.max(steady_h-272.0))
cbar = plt.colorbar(p, orientation='horizontal', ticks=np.arange(0.0,np.max(steady_h-272.0),5.0))
cbar.ax.set_xlabel('FLUID TEMPERATURE [$^{\circ}$C]')

# formatting
plt.yticks([0.0, -500.0, -1000.0], [0, -500, -1000.0])
plt.xticks([0.0, 1500.0, 3000.0], [0, 1500, 3000])
plt.xlim(np.min(x), np.max(x))


plt.savefig(path + 'plot_circulation.png')

    
for i in range(0,steps,1): 
    #i=1


    #######################
    # FORMAT GRIDDED DATA #
    #######################

##    h = h0[i*len(y)-i:((i)*len(y)+len(x))-i-1,:]
##    h = np.append(h, h[-1:,:], axis=0)
##    h = np.append(h, h[:,-1:], axis=1)
##
##    psi = psi0[i*len(y)-i:((i)*len(y)+len(x))-i-1,:]
##    psi = np.append(psi, psi[-1:,:], axis=0)
##    psi = np.append(psi, psi[:,-1:], axis=1)
##
##    perm = np.append(perm0, perm0[-1:,:], axis=0)
##    perm = np.append(perm, perm[:,-1:], axis=1)
##
##    v = v0[i*len(y)-i:((i)*len(y)+len(x))-i-1,:]
##    v = np.append(v, v[-1:,:], axis=0)
##    v = np.append(v, v[:,-1:], axis=1)
##
##    u = u0[i*len(y)-i:((i)*len(y)+len(x))-i-1,:]
##    u = np.append(u, u[-1:,:], axis=0)
##    u = np.append(u, u[:,-1:], axis=1)

    feldspar = feldspar0[(i*len(y0)/cell):(i*len(y0)/cell+len(y0)/cell)-1,:]
    feldspar = np.append(feldspar, feldspar[-1:,:], axis=0)
    feldspar = np.append(feldspar, feldspar[:,-1:], axis=1)

    glass = glass0[(i*len(y0)/cell):(i*len(y0)/cell+len(y0)/cell)-1,:]
    glass = np.append(glass, glass[-1:,:], axis=0)
    glass = np.append(glass, glass[:,-1:], axis=1)

    geo = geo0[(i*len(y0)/cell):(i*len(y0)/cell+len(y0)/cell),:]
    geo = np.append(geo, geo[-1:,:], axis=0)
    geo = np.append(geo, geo[:,-1:], axis=1)


    ###############
    # SECOND PLOT #
    ###############

    fig=plt.figure()
    ax1=fig.add_subplot(1,1,1, aspect='equal')

    #contours = np.round(np.arange(0.0,np.max(ca0),np.max(ca0)/10.0),6)

    
    #geoContours = np.arange(np.min(geo0),np.max(geo00)+(np.max(geo00)-np.min(geo0))/20.0,
                    #(np.max(geo00)-np.min(geo0))/20.0)
    #geoContours = np.linspace(0.00,0.0023,20)

    
    #ticks=np.arange(0.000,0.0028,.0004)
    #contours = np.arange(0.0020,0.0032,0.0004)
    #contours = np.arange(0.0,.162,0.02)
    # np.arange(0.0,0.0032,0.0004),
    print xCell.shape
    print yCell.shape
    print geo.shape
    #print yCell[:-1].shape

    # plot chem
    pGlass = plt.pcolor(xCell, yCell, geo, cmap=cm.autumn, linewidth=1.0, color='black')
    contoursPsi = np.arange(np.min(steady_psi),np.max(steady_psi)+(np.max(steady_psi)-np.min(steady_psi))/10.0,
                         (np.max(steady_psi)-np.min(steady_psi))/10.0)
    CS = plt.contour(xg, yg, steady_psi, contoursPsi, colors='#003399',linewidths=np.array([1.5]))


    # formatting
    #theTicks = geoContours
    cbar= plt.colorbar(pGlass, orientation='horizontal')
    cbar.ax.set_xlabel('DISSOLVED INORGANIC CARBON CONCENTRATION [mol/kgw]')
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    
    plt.title('t = ' + str(i*160) + ' years')

    
    plt.savefig(path + 'plot_e4_' + str(i) + '.png')




print "ALL DONE!"

