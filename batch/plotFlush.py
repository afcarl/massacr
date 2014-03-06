# plotFlush.py

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
from scipy.optimize import curve_fit

infile = 'r90t06c350pw.txt'
flush = np.loadtxt(infile)

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
# 13 hco3-
# 14 co32-
######### nothing al (NOW ALKALINITY)
# 15 stilbite
# 16
# 17 sio2am
# 18
# 19 kaolinite
# 20
# 21 albite
# 22
# 23 saponite-mg
# 24
# 25 celadonite
# 26
# 27 clinop
# 28
# 29 pyrite
# 30
# 31 hematite
# 32
# 33 goethite
# 34
# 35 dolomite
# 36
# 37 smectite
# 38
# 39 dawsonite
# 40
# 41 magnesite
# 42
# 43 siderite
# 44
# 45 calcite
# 46
# 47 quartz
# 48
# 49 k-spar
# 50
# 51 plagioclase
# 52
# 53 augite
# 54
# 55 pigeonite
# 55
# 57 magnetite
# 58
# 59 basaltic glass
# 60
# 61 phi
# 62 ssp
# 63 water volume
# 64 rho_s


print flush.shape

##print np.max(flush[:,51])-np.min(flush[:,51]), 'plag'
##print np.max(flush[:,53])-np.min(flush[:,53]), 'aug'
##print np.max(flush[:,55])-np.min(flush[:,55]), 'pig'
##print np.max(flush[:,57])-np.min(flush[:,57]), 'mag'
print flush[:,45], 'calcite'
print flush[:,13]
print flush[:,14]


flush[:,0] = flush[:,0]*flush[0,1]/(3.14e7)


fig=plt.figure()

fig.suptitle(infile, fontsize=14, fontweight='bold')

plt.rc('xtick', labelsize=6) 
plt.rc('ytick', labelsize=6)


#####################
# PLOT  DISSOLUTION #
#####################

ax = plt.subplot(2,2,1)
p = plt.plot(flush[:,0],flush[:,51],'r', label="plagioclase")
p = plt.plot(flush[:,0],flush[:,53],'g', label="augite")
p = plt.plot(flush[:,0],flush[:,55],'m--', label="pigeonite")
p = plt.plot(flush[:,0],flush[:,57],'gold', label="magnetite")
p = plt.plot(flush[:,0],flush[:,59],'k', label="basaltic glass")



handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.legend(handles, labels,loc='best',prop={'size':5}, ncol=2)

plt.ylabel('amount [mol]',fontsize=6)
plt.xlabel('time [yrs]',fontsize=6)


######################
# PLOT  PRECIPITATES #
######################

ax = plt.subplot(2,2,2)
p = plt.plot(flush[:,0],flush[:,15],'r', label="stilbite")
p = plt.plot(flush[:,0],flush[:,17],'g', label="sio2")
p = plt.plot(flush[:,0],flush[:,19],'b', label="kaolinite")
p = plt.plot(flush[:,0],flush[:,21],'m', label="albite")
p = plt.plot(flush[:,0],flush[:,23],'gold', label="saponite")
p = plt.plot(flush[:,0],flush[:,25],'grey', label="celadonite")
p = plt.plot(flush[:,0],flush[:,37],'k--', label="smectite")
p = plt.plot(flush[:,0],flush[:,27],'r--', label="Clinoptilolite-Ca")
p = plt.plot(flush[:,0],flush[:,49],'b--', label="k-spar")
#p = plt.plot(flush[:,0],flush[:,47],'purple', label="quartz")
p = plt.plot(flush[:,0],flush[:,45],'c-', label="calcite")

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
