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


mass = out[:,35]*270.0 + out[:,37]*230.0 + out[:,39]*238.6 +  \
       out[:,41]*231.0 + out[:,43]*46.5 + \
       out[:,3]*832.2 + out[:,5]*60.0 + out[:,7]*258.2 + out[:,9]*262.3 + \
        + out[:,11]*358.537 + out[:,13]*396.8 + \
       out[:,15]*943.16 + out[:,17]*120.0 + out[:,19]*159.6 + out[:,33]*130.0
#mass = 1.0


###############
# dissolution #
###############

ax = plt.subplot(2,2,2)



out[:,0] = out[:,0]/(3.14e7)
p3, = ax.plot(out[:,0],out[:,35]*270.0/mass, 'r', label='Plagioclase',linewidth=2)
p3, = ax.plot(out[:,0],out[:,37]*230.0/mass, 'g', label='Augite',linewidth=2)
p3, = ax.plot(out[:,0],out[:,39]*238.6/mass, 'm--', label='Pigeonite',linewidth=2)
p3, = ax.plot(out[:,0],out[:,41]*231.0/mass, 'gold', label='Magnetite',linewidth=2)
p3, = ax.plot(out[:,0],out[:,43]*46.5/mass, 'k', label='Basaltic Glass',linewidth=2)
plt.grid(True)



handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.yticks(np.arange(0,.5,.05))
plt.ylabel('MASS FRACTION',fontsize=8)
plt.xlabel('TIME',fontsize=8)
plt.legend(handles, labels,loc='best',prop={'size':8})

###########
# pH plot #
###########

ax = plt.subplot(2,2,1)
this = out[:,1]
p3, = ax.plot(out[:,0],out[:,1], 'r', label='pH',linewidth=2)
plt.grid(True)

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])

plt.ylim(3,8)
plt.yticks(np.arange(3,8,.5))
plt.ylabel('pH',fontsize=8)
plt.xlabel('TIME',fontsize=8)
plt.legend(handles, labels,loc='best',prop={'size':8})


################
# precipitates #
################

ax = plt.subplot(2,2,3)
this = out[:,1]
p3, = ax.plot(out[:,0],out[:,3]*832.2/mass, 'b:', label='Stilbite',linewidth=1)
p3, = ax.plot(out[:,0],out[:,5]*60.0/mass, 'b', label='SiO2',linewidth=1)
p3, = ax.plot(out[:,0],out[:,7]*258.2/mass, 'k', label='Kaolinite',linewidth=1)
p3, = ax.plot(out[:,0],out[:,9]*262.3/mass, 'g', label='Albite',linewidth=1)
p3, = ax.plot(out[:,0],out[:,11]*385.537/mass, 'm--', label='Saponite-Mg',linewidth=1)
p3, = ax.plot(out[:,0],out[:,13]*396.8/mass, 'r', label='Celadonite',linewidth=1)
p3, = ax.plot(out[:,0],out[:,15]*943.16/mass, 'm', label='Clinoptilolite-Ca',linewidth=1)
p3, = ax.plot(out[:,0],out[:,17]*120.0/mass, 'r--', label='Pyrite',linewidth=1)
p3, = ax.plot(out[:,0],out[:,19]*159.6/mass, 'gold', label='Hematite',linewidth=1)
##p3, = ax.plot(out[:,0],out[:,21]*230.0/mass, 'k:', label='',linewidth=1)
##p3, = ax.plot(out[:,0],out[:,23]*230.0/mass, 'k:', label='',linewidth=1)
##p3, = ax.plot(out[:,0],out[:,25]*230.0/mass, 'k:', label='',linewidth=1)
##p3, = ax.plot(out[:,0],out[:,27]*230.0/mass, 'k:', label='',linewidth=1)
##p3, = ax.plot(out[:,0],out[:,29]*230.0/mass, 'k:', label='',linewidth=1)
##p3, = ax.plot(out[:,0],out[:,31]*230.0/mass, 'k:', label='',linewidth=1)


handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])

plt.ylabel('MASS FRACTION',fontsize=8)
plt.xlabel('TIME',fontsize=8)
plt.legend(handles, labels,loc='best',prop={'size':6}, ncol=2)


################
#  carbonates  #
################

ax = plt.subplot(2,2,4)
this = out[:,1]
p3, = ax.plot(out[:,0],out[:,33]*150.0/mass, 'm--', label='Siderite',linewidth=1)

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])

plt.ylabel('MASS FRACTION',fontsize=8)
plt.xlabel('TIME',fontsize=8)
plt.legend(handles, labels,loc='best',prop={'size':8})


plt.savefig('basalt40c.eps')
