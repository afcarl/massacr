#plotGlass.py

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
from scipy.optimize import curve_fit


outPlag = np.loadtxt('outPlagioclase.txt')
outAug = np.loadtxt('outAugite.txt')
outPig = np.loadtxt('outPigeonite.txt')

out = np.loadtxt('outmat.txt')

def func(x, a, b, c, d, e):
    return a + b*x + c/x + d*np.log10(x) + e/(x*x)
x = outPlag[:,0]
yn = outPlag[:,1]
popt, pcov = curve_fit(func, x, yn)
print popt
yn = outAug[:,1]
popt, pcov = curve_fit(func, x, yn)
print popt
yn = outPig[:,1]
popt, pcov = curve_fit(func, x, yn)
print popt

fig=plt.figure()

plt.rc('xtick', labelsize=12) 
plt.rc('ytick', labelsize=12) 


ax = plt.subplot(1,1,1)
out[:,0] = np.round(out[:,0]/(3.14e7),2)
p3, = ax.plot(out[:,0],out[:,1], 'grey', label='plagioclase',linewidth=2)
p3, = ax.plot(out[:,0],out[:,3], 'r', label='augite',linewidth=2)
p3, = ax.plot(out[:,0],out[:,5], 'g', label='pigeonite',linewidth=2)
p3, = ax.plot(out[:,0],out[:,7], 'b', label='magnetite',linewidth=2)


handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.margins(0.05)
plt.xlabel('AMOUNT [mol]',fontsize=12)
plt.ylabel('TIME STEP',fontsize=12)
plt.legend(handles, labels,loc=1,prop={'size':12})




plt.savefig('glass.png')
