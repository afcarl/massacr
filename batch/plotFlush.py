# plotFlush.py

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
from scipy.optimize import curve_fit

flush = np.loadtxt('flush.txt')

# 0 step #
# 1 time (2)
# 2 ph (3)
# 3 ca (4)
# 4 mg
# 5 na
# 6 k
# 7 fe
# 8 s6
# 9 si
# 10 cl
# 11 al
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



print flush[:,60]




fig=plt.figure()

plt.rc('xtick', labelsize=6) 
plt.rc('ytick', labelsize=6)


##############
# PLOT  #
##############

ax = plt.subplot(2,2,1)
p = plt.plot(flush[:,0],flush[:,56],'o--')

plt.ylabel(' [mol]',fontsize=6)
plt.xlabel('time',fontsize=6)





plt.subplots_adjust(hspace=.4, wspace=.4)
plt.savefig('flush.png')
