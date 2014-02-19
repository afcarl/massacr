# plotFlush.py

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
from scipy.optimize import curve_fit

flush = np.loadtxt('continuous.txt')

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
######### nothing al (NOW ALKALINITY)
# 12 stilbite
# 13
# 14 sio2am
# 15
# 16 kaolinite
# 17
# 18 albite
# 19
# 20 saponite-mg
# 21
# 22 celadonite
# 23
# 24 clinop
# 25
# 26 pyrite
# 27
# 28 hematite
# 29
# 30 goethite
# 31
# 32 dolomite
# 33
# 34 smectite
# 35
# 36 dawsonite
# 37
# 38 magnesite
# 39
# 40 siderite
# 41
# 42 calcite
# 43
# 44 quartz
# 45
# 46 k-spar
# 47
# 48 plagioclase
# 49
# 50 augite
# 51
# 52 pigeonite
# 53
# 54 magnetite
# 55
# 56 basaltic glass
# 57
# 58 phi
# 59 ssp
# 60 water volume
# 61 rho_s


# to only look at timestep right before a flush
print flush.shape
#flush = flush[range(9,100,10),:]

flush[:,0] = flush[:,0]*flush[0,1]/(3.14e7)

print flush[:,0]
print flush[:,42]
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
p = plt.plot(flush[:,0],flush[:,48],'r', label="plagioclase")
p = plt.plot(flush[:,0],flush[:,50],'g', label="augite")
p = plt.plot(flush[:,0],flush[:,52],'m--', label="pigeonite")
p = plt.plot(flush[:,0],flush[:,54],'gold', label="magnetite")
p = plt.plot(flush[:,0],flush[:,56],'k', label="basaltic glass")

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.legend(handles, labels,loc='best',prop={'size':5}, ncol=2)

plt.ylabel('amount [mol]',fontsize=6)
plt.xlabel('time [yrs]',fontsize=6)


######################
# PLOT  PRECIPITATES #
######################

ax = plt.subplot(2,2,2)
p = plt.plot(flush[:,0],flush[:,12],'r', label="stilbite")
p = plt.plot(flush[:,0],flush[:,14],'g', label="sio2")
p = plt.plot(flush[:,0],flush[:,16],'b', label="kaolinite")
p = plt.plot(flush[:,0],flush[:,18],'m', label="albite")
p = plt.plot(flush[:,0],flush[:,20],'gold', label="saponite")
p = plt.plot(flush[:,0],flush[:,22],'grey', label="celadonite")
#p = plt.plot(flush[:,0],flush[:,44],'purple', label="quartz")
p = plt.plot(flush[:,0],flush[:,34],'k--', label="smectite")
p = plt.plot(flush[:,0],flush[:,24],'r--', label="Clinoptilolite-Ca")
p = plt.plot(flush[:,0],flush[:,46],'b--', label="k-spar")
#p = plt.plot(flush[:,0],flush[:,44],'purple', label="quartz")

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.legend(handles, labels,loc='best',prop={'size':5}, ncol=2)

plt.ylabel('amount [mol]',fontsize=6)
plt.xlabel('time [yrs]',fontsize=6)



####################
# PLOT  CARBONATES #
####################

ax = plt.subplot(2,2,3)
#p = plt.plot(flush[:,0],flush[:,44],'r', label="quartz")
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
