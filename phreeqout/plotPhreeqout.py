#plotPhreeqout.py

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math

outmat = np.loadtxt('outmat.txt')

x = outmat[:,0]
y = outmat[:,1]

fig=plt.figure()

ax = plt.subplot(1,1,1)
#p1, = ax.plot(outmat[:,0],outmat[:,1], label='QUARTZ')
p2, = ax.plot(outmat[:,0],outmat[:,1], label='Ca$^{2+}$',linewidth=4)
p3, = ax.plot(outmat[:,0],outmat[:,2], label='Mg$^{2+}$',linewidth=4)
#p4, = ax.plot(outmat[:,0],outmat[:,5], label='DOLOMITE')

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
#plt.xlim(np.min(x)-abs(x[0]-x[-1])/10.0,np.max(x)+abs(x[0]-x[-1])/10.0)
#plt.ylim(np.min(y)-abs(y[0]-y[-1])/10.0,0.05)
plt.margins(0.05)
plt.xlabel('TEMPERATURE [C]')
plt.ylabel('MOLALITY IN SEAWATER [mol/kg]')
plt.legend(handles, labels,loc=3)
plt.savefig('ga.png')
