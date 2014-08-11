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

cell = 4
steps = 400
#path = "output/noTransportCell2/"
#path = "output/noTransportCell2/"
path = ""
#path = "output/transportWorks1200yr/"
#path = "output/par4cell1/"
#path = "output/par_4_2_cell1/"
#path = "output/sixty/"
#path = "output/off_test/"
#path = "output/test6transport/"
#path = "output/testLR/"
#path = "output/interpBlock/"
#path = "output/fix/"
#path = "output/perm/"

t = np.loadtxt(path + 't.txt',delimiter='\n')
x0 = np.loadtxt(path + 'x.txt',delimiter='\n')
y0 = np.loadtxt(path + 'y.txt',delimiter='\n')

x=x0
y=y0
bits = len(x)
x = np.append(x0, np.max(x0)+.001)
y = np.append(y0, np.max(y0)+(y0[-1]-y0[-2]))


xCell = x0
yCell = y0
xCell = xCell[::cell]
yCell= yCell[::cell]

xCell = np.append(xCell, np.max(xCell)+xCell[0])
yCell = np.append(yCell, np.max(yCell)-yCell[-1])

xg, yg = np.meshgrid(x[:],y[:])

h0 = np.loadtxt(path + 'hMat.txt')
u0= np.loadtxt(path + 'uMat.txt')
v0= np.loadtxt(path + 'vMat.txt')
psi0 = np.loadtxt(path + 'psiMat.txt')
feldspar0 = np.loadtxt(path + 'pri_feldspar.txt') 
glass0 = np.loadtxt(path + 'pri_glass.txt')
perm0 = np.loadtxt(path + 'permeability.txt')
geo0 = np.loadtxt(path + 'sol_c.txt')


geo00 = np.zeros(steps)
print "go"
for i in range(steps):
    geo1 = geo0[(i*len(y0)/cell):(i*len(y0)/cell+len(y0)/cell)-1,:]
    geo1 = np.append(geo1, geo1[-1:,:], axis=0)
    geo1 = np.append(geo1, geo1[:,-1:], axis=1)
    geo00[i] = np.max(geo1)
    
for i in range(0,steps,5): 
    #i=1
    print h0.shape


    #######################
    # FORMAT GRIDDED DATA #
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

    geo = geo0[(i*len(y0)/cell):(i*len(y0)/cell+len(y0)/cell)-1,:]
    geo = np.append(geo, geo[-1:,:], axis=0)
    geo = np.append(geo, geo[:,-1:], axis=1)


    ##############
    # FIRST PLOT #
    ##############

    fig=plt.figure()

    ax1=fig.add_subplot(1,1,1, aspect='equal')
    #levels00 = np.linspace(.000002, np.max(psi), 15)
    #levels0 = np.linspace(np.min(psi), -.000002, 15)
    #levels = np.append(levels0,levels00,axis=1)

    # permeability plot
    permC = plt.contour(xg, yg, np.log10(perm), [-14.0,-14.1], colors='w',linewidths=np.array([2.0]))
    #permC = plt.contourf(xg, yg, np.log10(perm), 10, cmap=cm.summer)

    # levels[::2],
    #plt.clabel(CS,  inline=0, fmt='>', fontsize=14)
    CS = plt.contour(xg, yg, psi, 10, colors='k',linewidths=np.array([1.0]))

    p = plt.contourf(xg,yg,h-272.0, np.arange(0.0,np.max(h-272.0),5.0), cmap=cm.rainbow)
    plt.clim(0.0,np.max(h-272.0))
    cbar = plt.colorbar(p, orientation='horizontal', ticks=np.arange(0.0,np.max(h-272.0),5.0))
    cbar.ax.set_xlabel('FLUID TEMPERATURE [$^{\circ}$C]')

    plt.yticks([0.0, -500.0, -1000.0], [0, -500, -1000.0])
    plt.xticks([0.0, 1500.0, 3000.0], [0, 1500, 3000])

    #plt.title("STREAMFUNCTIONS",fontsize=8)

    plt.xlim(np.min(x), np.max(x))

    plt.savefig('j14.png')


    ###############
    # SECOND PLOT #
    ###############

    fig=plt.figure()

    ax1=fig.add_subplot(1,1,1, aspect='equal')

    #contours = np.round(np.arange(0.0,np.max(ca0),np.max(ca0)/10.0),6)
    contours = np.arange(np.min(geo0),np.max(geo00)+(np.max(geo00)-np.min(geo0))/10.0,
                         (np.max(geo00)-np.min(geo0))/10.0)
    #ticks=np.arange(0.000,0.0028,.0004)
    contours = np.arange(0.0,0.0032,0.0004)
    #contours = np.arange(6.0,9.0,0.25)
    print geo.shape
    print xCell.shape
    print yCell[:-1].shape
    # np.arange(0.0,0.0032,0.0004),
    pGlass = plt.contourf(xCell, yCell[:-1], geo, np.arange(0.0,0.0032,0.0004), cmap=cm.YlOrRd)

    #FF6600
    contoursPsi = np.arange(np.min(psi),np.max(psi)+(np.max(psi)-np.min(psi))/10.0,
                         (np.max(psi)-np.min(psi))/10.0)
    CS = plt.contour(xg, yg, psi, contoursPsi, colors='#003399',linewidths=np.array([1.5]))

    plt.clabel(CS, contoursPsi, inline=0, fmt='  >  ', fontsize=14,
               rightside_up="True")

    theTicks = contours
    cbar= plt.colorbar(pGlass, orientation='horizontal')
    cbar.ax.set_xlabel('DISSOLVED INORGANIC CARBON [mol/kgw]')
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    #ticks=np.arange(0.0,0.0045,0.0009)
    #cbar.set_clim(vmin=0.0,vmax=.012)

    plt.title('t = ' + str(i*128) + ' years')

    plt.savefig(path + 'cm0' + str(i) + '.png')


print "ALL DONE!"

