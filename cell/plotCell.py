# plotCell.py

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
from scipy.optimize import curve_fit

infile = 'cell3.txt'
flush = np.loadtxt(infile)


# 0 step #
# 1 time (2)
# 2 ph (3)
# 3 pe (4)
# 4 alkalinity (5)
# 5 total C (6)
# 6 ca 
# 7 mg
# 8 na
# 9 k
# 10 fe
# 11 s6
# 12 si
# 13 cl
# 14 al
# 15 hco3-
###### 14 co32-
######### nothing al (NOW ALKALINITY)
# 16 stilbite (17)
# 17
# 18 sio2am
# 19
# 20 kaolinite
# 21
# 22 albite
# 23
# 24 saponite-mg
# 25
# 26 celadonite
# 27
# 28 clinop
# 29
# 30 pyrite
# 31
# 32 hematite / Montmor-Na
# 33
# 34 goethite
# 35
# 36 dolomite / boehmite
# 37
# 38 smectite
# 39
# 40 dawsonite
# 41
# 42 magnesite / analcime
# 43
# 44 siderite
# 45
# 46 calcite
# 47
# 48 quartz
# 49
# 50 k-spar
# 51
# 52 Saponite-Na
# 53
# 54 Nontronite-Na
# 55
# 56 Nontronite-Mg
# 57
# 58 Nontronite-K
# 59
# 60 Nontronite-H
# 61
# 62 Nontronite-Ca
# 63
# 64 muscovite
# 65
# 66 mesolite
# 67
# 68 hematite
# 69
# 70 diaspore (71)
# 71 
# 72 plagioclase (73)
# 73
# 74 augite
# 75
# 76 pigeonite
# 77
# 78 magnetite
# 79
# 80 basaltic glass (81)
# 81
# 82 phi
# 83 ssp
# 84 water volume
# 85 rho_s


print flush.shape

print flush[:,46], 'calcite'


flush[:,0] = flush[:,0]*flush[0,1]/(3.14e7)


fig=plt.figure()

fig.suptitle(infile, fontsize=14, fontweight='bold')

plt.rc('xtick', labelsize=6) 
plt.rc('ytick', labelsize=6)


#####################
# PLOT  DISSOLUTION #
#####################

ax = plt.subplot(2,2,1)
p = plt.plot(flush[:,0],flush[:,72],'r', label="plagioclase")
p = plt.plot(flush[:,0],flush[:,74],'g', label="augite")
p = plt.plot(flush[:,0],flush[:,76],'m--', label="pigeonite")
p = plt.plot(flush[:,0],flush[:,78],'gold', label="magnetite")
p = plt.plot(flush[:,0],flush[:,80],'k', label="basaltic glass")



handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.legend(handles, labels,loc='best',prop={'size':5}, ncol=2)

plt.ylabel('amount [mol]',fontsize=6)
plt.xlabel('time [yrs]',fontsize=6)


######################
# PLOT  PRECIPITATES #
######################

ax = plt.subplot(2,2,2)
p = plt.plot(flush[:,0],flush[:,16], label="stilbite")
p = plt.plot(flush[:,0],flush[:,18], label="sio2")
p = plt.plot(flush[:,0],flush[:,20], label="kaolinite")
p = plt.plot(flush[:,0],flush[:,22], label="albite")
p = plt.plot(flush[:,0],flush[:,24], label="saponite")
p = plt.plot(flush[:,0],flush[:,26], label="celadonite")
p = plt.plot(flush[:,0],flush[:,28], label="clinop")
#p = plt.plot(flush[:,0],flush[:,30], label="pyrite")
p = plt.plot(flush[:,0],flush[:,32], '--', label="mont-na")
p = plt.plot(flush[:,0],flush[:,34], '--', label="goethite")
p = plt.plot(flush[:,0],flush[:,36], '--', label="dolomite")
p = plt.plot(flush[:,0],flush[:,38], '--', label="smectite")
p = plt.plot(flush[:,0],flush[:,40], '--', label="dawsonite")
p = plt.plot(flush[:,0],flush[:,42], '--', label="magnesite")
#p = plt.plot(flush[:,0],flush[:,44], label="siderite")
p = plt.plot(flush[:,0],flush[:,46], '--', label="calcite")
#p = plt.plot(flush[:,0],flush[:,48], ':', label="qtz")
p = plt.plot(flush[:,0],flush[:,50], ':', label="k-spar")
p = plt.plot(flush[:,0],flush[:,52], ':', label="saponite-na")
p = plt.plot(flush[:,0],flush[:,54], ':', label="nont-na")
p = plt.plot(flush[:,0],flush[:,56], ':', label="nont-mg")
p = plt.plot(flush[:,0],flush[:,58], ':', label="nont-k")
p = plt.plot(flush[:,0],flush[:,60], ':', label="nont-h")
p = plt.plot(flush[:,0],flush[:,62], linewidth=2, label="nont-ca")
p = plt.plot(flush[:,0],flush[:,64], linewidth=2, label="muscovite")
p = plt.plot(flush[:,0],flush[:,66], linewidth=2, label="mesolite")
p = plt.plot(flush[:,0],flush[:,68], linewidth=2, label="hematite")
p = plt.plot(flush[:,0],flush[:,70], linewidth=2, label="diaspore")

##p = plt.plot(flush[:,0],flush[:,38], label="smectite")
##p = plt.plot(flush[:,0],flush[:,28], label="Clinoptilolite-Ca")
##p = plt.plot(flush[:,0],flush[:,50], label="k-spar")
##p = plt.plot(flush[:,0],flush[:,32], label="Montmor-Na",linewidth=2)
##p = plt.plot(flush[:,0],flush[:,46], label="calcite")
##p = plt.plot(flush[:,0],flush[:,36], label="dolomite")

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.legend(handles, labels,loc='best',prop={'size':5}, ncol=2)

plt.ylabel('amount [mol]',fontsize=6)
plt.xlabel('time [yrs]',fontsize=6)




# pH #


ax = plt.subplot(2,2,3)
#p = plt.plot(flush[:,0],flush[:,46],'r', label="quartz")
p = plt.plot(flush[:,0],flush[:,2],'g', label="ph")

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.legend(handles, labels,loc='best',prop={'size':5}, ncol=2)


plt.ylabel('amount [mol]',fontsize=6)
plt.xlabel('time [yrs]',fontsize=6)





# ALKALINITY #


ax = plt.subplot(2,2,4)
p = plt.plot(flush[:,0],flush[:,4],'g', label="alk")

#+2.0*flush[:,35]/flush[:,83]

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.legend(handles, labels,loc='best',prop={'size':5}, ncol=2)


plt.ylabel('amount [mol]',fontsize=6)
plt.xlabel('time [yrs]',fontsize=6)



plt.subplots_adjust(hspace=.25, wspace=.25)
plt.savefig('cell.png')
