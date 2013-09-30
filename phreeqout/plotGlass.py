#plotGlass.py

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
from scipy.optimize import curve_fit


outPlag = np.loadtxt('outPlagioclase.txt')
outAug = np.loadtxt('outAugite.txt')
outPig = np.loadtxt('outPigeonite.txt')

out = np.loadtxt('testMat.txt')

def func(x, a, b, c, d, e):
    return a + b*x + c/x + d*np.log10(x) + e/(x*x)
x = outPlag[:,0]+273.0
yn = outPlag[:,1]
popt, pcov = curve_fit(func, x, yn)
print popt
yn = outAug[:,1]
popt, pcov = curve_fit(func, x, yn)
print popt
yn = outPig[:,1]
popt, pcov = curve_fit(func, x, yn)
print popt

print outPlag[11,1]
print outAug[11,1]
print outPig[11,1]

fig=plt.figure()

plt.rc('xtick', labelsize=8) 
plt.rc('ytick', labelsize=8)


mass = out[:,34]*270.0 + out[:,36]*230.0 + out[:,38]*238.6 +  \
       out[:,40]*231.0 + out[:,42]*46.5 + \
       out[:,2]*817.2 + out[:,4]*60.0 + out[:,6]*258.2 + out[:,8]*262.3 + \
        + out[:,10]*358.537 + out[:,12]*396.8 + \
       out[:,14]*943.16 + out[:,16]*120.0 + out[:,18]*159.6
#mass = 1.0
print mass

###############
# dissolution #
###############

ax = plt.subplot(2,2,2)



out[:,0] = np.round(out[:,0]/(3.14e7),2)
p3, = ax.plot(out[:,0],out[:,34]*270.0/mass, 'r', label='plagioclase',linewidth=2)
p3, = ax.plot(out[:,0],out[:,36]*230.0/mass, 'g', label='augite',linewidth=2)
p3, = ax.plot(out[:,0],out[:,38]*238.6/mass, 'm--', label='pigeonite',linewidth=2)
p3, = ax.plot(out[:,0],out[:,40]*231.0/mass, 'gold', label='magnetite',linewidth=2)
p3, = ax.plot(out[:,0],out[:,42]*46.5/mass, 'k', label='basaltic glass',linewidth=2)
plt.grid(True)

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
#plt.margins(0.05)
plt.ylabel('MASS FRACTION',fontsize=8)
plt.xlabel('TIME ELAPSED [year]',fontsize=8)
plt.legend(handles, labels,loc='best',prop={'size':8})

###########
# pH plot #
###########

ax = plt.subplot(2,2,1)
this = out[:,1]
p3, = ax.plot(out[:,0],out[:,1], 'grey', label='pH',linewidth=2)

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.margins(0.05)
plt.ylabel('pH',fontsize=8)
plt.xlabel('TIME [year]',fontsize=8)
plt.legend(handles, labels,loc='best',prop={'size':8})


################
# precipitates #
################

ax = plt.subplot(2,2,4)
this = out[:,1]
p3, = ax.plot(out[:,0],out[:,2]*817.2/mass, 'b:', label='stilbite',linewidth=1)
p3, = ax.plot(out[:,0],out[:,4]*60.0/mass, 'b', label='sio2',linewidth=1)
p3, = ax.plot(out[:,0],out[:,6]*258.2/mass, 'k', label='kaolinite',linewidth=1)
p3, = ax.plot(out[:,0],out[:,8]*262.3/mass, 'g', label='albite',linewidth=1)
p3, = ax.plot(out[:,0],out[:,10]*385.537/mass, 'm--', label='saponite-mg',linewidth=1)
p3, = ax.plot(out[:,0],out[:,12]*396.8/mass, 'r', label='celadonite',linewidth=1)
p3, = ax.plot(out[:,0],out[:,14]*943.16/mass, 'm', label='Clinoptilolite-Ca',linewidth=1)
p3, = ax.plot(out[:,0],out[:,16]*120.0/mass, 'r--', label='pyrite',linewidth=1)
p3, = ax.plot(out[:,0],out[:,18]*159.6/mass, 'gold', label='hematite',linewidth=1)
##p3, = ax.plot(out[:,0],out[:,20]*230.0/mass, 'k:', label='',linewidth=1)
##p3, = ax.plot(out[:,0],out[:,22]*230.0/mass, 'k:', label='',linewidth=1)
##p3, = ax.plot(out[:,0],out[:,24]*230.0/mass, 'k:', label='',linewidth=1)
##p3, = ax.plot(out[:,0],out[:,26]*230.0/mass, 'k:', label='',linewidth=1)
##p3, = ax.plot(out[:,0],out[:,28]*230.0/mass, 'k:', label='',linewidth=1)
##p3, = ax.plot(out[:,0],out[:,30]*230.0/mass, 'k:', label='',linewidth=1)


handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.margins(0.05)
plt.ylabel('MASS FRACTION',fontsize=8)
plt.xlabel('TIME [year]',fontsize=8)
plt.legend(handles, labels,loc='best',prop={'size':7})


################
# mass ? #
################

ax = plt.subplot(2,2,3)
this = out[:,1]
p3, = ax.plot(out[:,0],out[:,18]*159.6/mass, 'm--', label='siderite',linewidth=1)

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.margins(0.05)
plt.ylabel('moles total',fontsize=8)
plt.xlabel('TIME [year]',fontsize=8)
plt.legend(handles, labels,loc='best',prop={'size':8})


plt.savefig('glass927.png')
