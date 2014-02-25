# plotFlush.py

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
from scipy.optimize import curve_fit

flush = np.loadtxt('cmdTest.txt')

# 0 step #
# 1 time (2)
# 2 ph (3)
# 3 alkalinity (4)
# 4 ca 
# 5 mg
# 6 na
# 7 k
# 8 fe
# 9 s6
# 10 si
# 11 cl
# 12 al
######### nothing al (NOW ALKALINITY)
# 13 stilbite
# 14
# 15 sio2am
# 16
# 17 kaolinite
# 18
# 19 albite
# 20
# 21 saponite-mg
# 22
# 23 celadonite
# 24
# 25 clinop
# 26
# 27 pyrite
# 28
# 29 hematite
# 30
# 31 goethite
# 32
# 33 dolomite
# 34
# 35 smectite
# 36
# 37 dawsonite
# 38
# 39 magnesite
# 40
# 41 siderite
# 42
# 43 calcite
# 44
# 45 quartz
# 46
# 47 k-spar
# 48
# 49 plagioclase
# 50
# 51 augite
# 52
# 53 pigeonite
# 54
# 55 magnetite
# 56
# 57 basaltic glass
# 58
# 59 phi
# 60 ssp
# 61 water volume
# 62 rho_s


print flush.shape
flush[:,0] = flush[:,0]*flush[0,1]/(3.14e7)

##print "clintop"
##print flush[:,24]
##print "molal"
##print flush[:,5]
#print flush[:,14]
#print flush[:,16]
#print flush[:,18]


fig=plt.figure()

plt.rc('xtick', labelsize=6) 
plt.rc('ytick', labelsize=6)


#####################
# PLOT  DISSOLUTION #
#####################

ax = plt.subplot(2,2,1)
p = plt.plot(flush[:,0],flush[:,49],'r', label="plagioclase")
p = plt.plot(flush[:,0],flush[:,51],'g', label="augite")
p = plt.plot(flush[:,0],flush[:,53],'m--', label="pigeonite")
p = plt.plot(flush[:,0],flush[:,55],'gold', label="magnetite")
p = plt.plot(flush[:,0],flush[:,57],'k', label="basaltic glass")

print np.max(flush[:,49])-np.min(flush[:,49]), 'plag'
print np.max(flush[:,51])-np.min(flush[:,51]), 'aug'
print np.max(flush[:,53])-np.min(flush[:,53]), 'pig'
print np.max(flush[:,55])-np.min(flush[:,55]), 'mag'
print flush[:,43], 'calcite'

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.legend(handles, labels,loc='best',prop={'size':5}, ncol=2)

plt.ylabel('amount [mol]',fontsize=6)
plt.xlabel('time [yrs]',fontsize=6)


######################
# PLOT  PRECIPITATES #
######################

ax = plt.subplot(2,2,2)
p = plt.plot(flush[:,0],flush[:,13],'r', label="stilbite")
p = plt.plot(flush[:,0],flush[:,15],'g', label="sio2")
p = plt.plot(flush[:,0],flush[:,17],'b', label="kaolinite")
p = plt.plot(flush[:,0],flush[:,19],'m', label="albite")
p = plt.plot(flush[:,0],flush[:,21],'gold', label="saponite")
p = plt.plot(flush[:,0],flush[:,23],'grey', label="celadonite")
#p = plt.plot(flush[:,0],flush[:,45],'purple', label="quartz")
p = plt.plot(flush[:,0],flush[:,35],'k--', label="smectite")
p = plt.plot(flush[:,0],flush[:,25],'r--', label="Clinoptilolite-Ca")
p = plt.plot(flush[:,0],flush[:,47],'b--', label="k-spar")
#p = plt.plot(flush[:,0],flush[:,45],'purple', label="quartz")
p = plt.plot(flush[:,0],flush[:,43],'c-', label="calcite")

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.legend(handles, labels,loc='best',prop={'size':5}, ncol=2)

plt.ylabel('amount [mol]',fontsize=6)
plt.xlabel('time [yrs]',fontsize=6)



####################
# PLOT  CARBONATES #
####################

ax = plt.subplot(2,2,3)
#p = plt.plot(flush[:,0],flush[:,45],'r', label="quartz")
p = plt.plot(flush[:,0],flush[:,2],'g', label="ph")

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.legend(handles, labels,loc='best',prop={'size':5}, ncol=2)


plt.ylabel('amount [mol]',fontsize=6)
plt.xlabel('time [yrs]',fontsize=6)




##############
# PLOT  IONS #
##############

ax = plt.subplot(2,2,4)
p = plt.plot(flush[:,0],flush[:,3],'g', label="alk")

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.legend(handles, labels,loc='best',prop={'size':5}, ncol=2)


plt.ylabel('amount [mol]',fontsize=6)
plt.xlabel('time [yrs]',fontsize=6)




plt.subplots_adjust(hspace=.25, wspace=.25)
plt.savefig('flush.png')
