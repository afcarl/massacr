#plotAlt.py

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math


out = np.loadtxt('outmat.txt')

# ODP LEG 168 DATA #

temp168 = [15.5,22.8,38.6,40.1,40.4,50.5,58.7,57.1,61.7,61.7,62.8]
mg168 = [47.23,42.51,27.33,5.66,4.65,16.02,12.67,2.79,4.64,2.21,4.98]
ca168 = [11.54,16.69,36.39,57.12,60.27,40.36,54.46,61.44,65.72,56.13,62.47]
alk168 = [2.89,1.83,0.92,0.46,0.39,0.38,0.79,0.27,0.35,0.63,0.48]

fig=plt.figure()

plt.rc('xtick', labelsize=10) 
plt.rc('ytick', labelsize=10) 

##################
# Ca / Mg switch #
##################


ax = plt.subplot(1,1,1)

p2, = ax.plot(out[:,0],1000.0*out[:,2], label='Ca$^{2+}$ (NAVAH)',linewidth=2)
p3, = ax.plot(out[:,0],1000.0*out[:,3], label='Mg$^{2+}$ (NAVAH)',linewidth=2)
p2, = ax.plot(temp168,ca168,'bo--', label='Ca$^{2+}$ (ODP LEG 168)',linewidth=2)
p3, = ax.plot(temp168,mg168,'go--', label='Mg$^{2+}$ (ODP LEG 168)',linewidth=2)

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.margins(0.05)
plt.xlabel('TEMPERATURE [C]',fontsize=12)
plt.ylabel('MOLALITY IN SEAWATER [mmol/kg]',fontsize=12)
plt.legend(handles, labels,loc=2,prop={'size':10})
plt.xlim(xmax=100)

##################
#  ALKALINITY    #
##################


##ax = plt.subplot(1,2,2)
##
##p2, = ax.plot(out[:,0],out[:,1], label='ALKALINITY (NAVAH',linewidth=2)
##p2, = ax.plot(temp168,alk168,'bo--', label='ALKALINITY (ODP LEG 168)',linewidth=2)
##
##
##handles, labels = ax.get_legend_handles_labels()
##plt.legend(handles[::-1], labels[::-1])
##plt.margins(0.05)
##plt.xlabel('TEMPERATURE [C]',fontsize=8)
##plt.ylabel('MOLALITY IN SEAWATER [mmol/kg]',fontsize=8)
##plt.legend(handles, labels,loc=1,prop={'size':6})


plt.savefig('leg168.png')
