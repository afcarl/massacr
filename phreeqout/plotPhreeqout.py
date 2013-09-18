#plotPhreeqout.py

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math




fig=plt.figure()

plt.rc('xtick', labelsize=6) 
plt.rc('ytick', labelsize=6) 

##################
# Ca / Mg switch #
##################

out = np.loadtxt('outForstDiop.txt')

ax = plt.subplot(2,3,1)

p2, = ax.plot(out[:,0],out[:,2], label='Ca$^{2+}$',linewidth=2)
p3, = ax.plot(out[:,0],out[:,3], label='Mg$^{2+}$',linewidth=2)


handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.margins(0.05)
#plt.xlabel('TEMPERATURE [C]',fontsize=8)
plt.ylabel('MOLALITY IN SEAWATER [mol/kg]',fontsize=8)
plt.legend(handles, labels,loc=1,prop={'size':6})
plt.title('Forsterite & Diopside',fontsize=8)


########################
#    EVERYTHING        #
########################

out = np.loadtxt('outForstWollast.txt')

ax = plt.subplot(2,3,2)

p2, = ax.plot(out[:,0],out[:,2], label='Ca$^{2+}$',linewidth=2)
p3, = ax.plot(out[:,0],out[:,3], label='Mg$^{2+}$',linewidth=2)


handles, labels = ax.get_legend_handles_labels()
#plt.legend(handles[::-1], labels[::-1])
plt.margins(0.05)
#plt.xlabel('TEMPERATURE [C]')
#plt.ylabel('MOLALITY IN SEAWATER [mol/kg]',fontsize=8)
#plt.legend(handles, labels,loc=1,prop={'size':8})
plt.title('Forsterite & Wollastonite',fontsize=8)

########################
#    EVERYTHING        #
########################

out = np.loadtxt('outForstWollastAnorth.txt')

ax = plt.subplot(2,3,3)

p2, = ax.plot(out[:,0],out[:,2], label='Ca$^{2+}$',linewidth=2)
p3, = ax.plot(out[:,0],out[:,3], label='Mg$^{2+}$',linewidth=2)


handles, labels = ax.get_legend_handles_labels()
#plt.legend(handles[::-1], labels[::-1])
plt.margins(0.05)
#plt.xlabel('TEMPERATURE [C]')
#plt.ylabel('MOLALITY IN SEAWATER [mol/kg]',fontsize=8)
#plt.legend(handles, labels,loc=1,prop={'size':8})
plt.title('Forsterite, Wollastonite & Anorthite',fontsize=8)



########################
#    EVERYTHING        #
########################

out = np.loadtxt('outQuartzCaol.txt')

ax = plt.subplot(2,3,4)

p2, = ax.plot(out[:,0],out[:,2], label='Ca$^{2+}$',linewidth=2)
p3, = ax.plot(out[:,0],out[:,3], label='Mg$^{2+}$',linewidth=2)


handles, labels = ax.get_legend_handles_labels()
#plt.legend(handles[::-1], labels[::-1])
plt.margins(0.05)
plt.xlabel('TEMPERATURE [C]',fontsize=8)
#plt.ylabel('MOLALITY IN SEAWATER [mol/kg]',fontsize=8)
#plt.legend(handles, labels,loc=1,prop={'size':8})
plt.title('Quartz & Ca-olivine',fontsize=8)

########################
#    EVERYTHING        #
########################

out = np.loadtxt('outQuartzForstDiop.txt')

ax = plt.subplot(2,3,5)

p2, = ax.plot(out[:,0],out[:,2], label='Ca$^{2+}$',linewidth=2)
p3, = ax.plot(out[:,0],out[:,3], label='Mg$^{2+}$',linewidth=2)


handles, labels = ax.get_legend_handles_labels()
#plt.legend(handles[::-1], labels[::-1])
plt.margins(0.05)
#plt.xlabel('TEMPERATURE [C]')
#plt.ylabel('MOLALITY IN SEAWATER [mol/kg]',fontsize=8)
#plt.legend(handles, labels,loc=1,prop={'size':8})
plt.title('Quartz, Forsterite & Diopside',fontsize=8)

########################
#    EVERYTHING        #
########################

out = np.loadtxt('outCaolWollast.txt')

ax = plt.subplot(2,3,6)

p2, = ax.plot(out[:,0],out[:,2], label='Ca$^{2+}$',linewidth=2)
p3, = ax.plot(out[:,0],out[:,3], label='Mg$^{2+}$',linewidth=2)


handles, labels = ax.get_legend_handles_labels()
#plt.legend(handles[::-1], labels[::-1])
plt.margins(0.05)
#plt.xlabel('TEMPERATURE [C]')
#plt.ylabel('MOLALITY IN SEAWATER [mol/kg]',fontsize=8)
#plt.legend(handles, labels,loc=1,prop={'size':8})
plt.title('Ca-olivine & Wollastonite',fontsize=8)


plt.savefig('assemblages.eps')
