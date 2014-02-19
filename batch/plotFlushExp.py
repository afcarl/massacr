# plotFlushExp.py

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
from scipy.optimize import curve_fit


temps = np.arange(2,22,2)


#############################
# SIMULATIONS WITH MINERALS #
# MIXING RATIO 90/10        #
#############################

# load temps
t02 = np.loadtxt('r90t02.txt')
t04 = np.loadtxt('r90t04.txt')
t06 = np.loadtxt('r90t06.txt')
t08 = np.loadtxt('r90t08.txt')
t10 = np.loadtxt('r90t10.txt')
t12 = np.loadtxt('r90t12.txt')
t14 = np.loadtxt('r90t14.txt')
t16 = np.loadtxt('r90t16.txt')
t18 = np.loadtxt('r90t18.txt')
t20 = np.loadtxt('r90t20.txt')

t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]]]

# grab for each experiment
t90_alk = np.zeros((10))
t90_glass = np.zeros((10))
t90_kaolinite = np.zeros((10))
t90_stilbite = np.zeros((10))
t90_saponite = np.zeros((10))
t90_albite = np.zeros((10))
t90_celadonite = np.zeros((10))
for i in range(len(temps)):
    bit = np.asarray(t[i])
    t90_alk[i] = bit[0,-2,3] - .002*.1
    t90_glass[i] = bit[0,-2,57]
    t90_kaolinite[i] = bit[0,-3,18]
    t90_stilbite[i] = bit[0,-3,14]
    t90_saponite[i] = bit[0,-3,22]
    t90_albite[i] = bit[0,-3,20]
    t90_celadonite[i] = bit[0,-3,24]


#############################
# SIMULATIONS WITH MINERALS #
# MIXING RATIO 80/20        #
#############################

# load temps
t02 = np.loadtxt('r80t02.txt')
t04 = np.loadtxt('r80t04.txt')
t06 = np.loadtxt('r80t06.txt')
t08 = np.loadtxt('r80t08.txt')
t10 = np.loadtxt('r80t10.txt')
t12 = np.loadtxt('r80t12.txt')
t14 = np.loadtxt('r80t14.txt')
t16 = np.loadtxt('r80t16.txt')
t18 = np.loadtxt('r80t18.txt')
t20 = np.loadtxt('r80t20.txt')

t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]]]

# grab for each experiment
t80_alk = np.zeros((10))
t80_glass = np.zeros((10))
t80_kaolinite = np.zeros((10))
t80_stilbite = np.zeros((10))
t80_saponite = np.zeros((10))
t80_albite = np.zeros((10))
t80_celadonite = np.zeros((10))
for i in range(len(temps)):
    bit = np.asarray(t[i])
    t80_alk[i] = bit[0,-2,3] - .002*.2
    t80_glass[i] = bit[0,-2,57]
    t80_kaolinite[i] = bit[0,-3,18]
    t80_stilbite[i] = bit[0,-3,14]
    t80_saponite[i] = bit[0,-3,22]
    t80_albite[i] = bit[0,-3,20]
    t80_celadonite[i] = bit[0,-3,24]





#############################
# SIMULATIONS WITH MINERALS #
# MIXING RATIO 70/30        #
#############################

# load temps
t02 = np.loadtxt('r70t02.txt')
t04 = np.loadtxt('r70t04.txt')
t06 = np.loadtxt('r70t06.txt')
t08 = np.loadtxt('r70t08.txt')
t10 = np.loadtxt('r70t10.txt')
t12 = np.loadtxt('r70t12.txt')
t14 = np.loadtxt('r70t14.txt')
t16 = np.loadtxt('r70t16.txt')
t18 = np.loadtxt('r70t18.txt')
t20 = np.loadtxt('r70t20.txt')

t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]]]

# grab for each experiment
t70_alk = np.zeros((10))
t70_glass = np.zeros((10))
t70_kaolinite = np.zeros((10))
t70_stilbite = np.zeros((10))
t70_saponite = np.zeros((10))
t70_albite = np.zeros((10))
t70_celadonite = np.zeros((10))
for i in range(len(temps)):
    bit = np.asarray(t[i])
    t70_alk[i] = bit[0,-2,3] - .002*.3
    t70_glass[i] = bit[0,-2,57]
    t70_kaolinite[i] = bit[0,-3,18]
    t70_stilbite[i] = bit[0,-3,14]
    t70_saponite[i] = bit[0,-3,22]
    t70_albite[i] = bit[0,-3,20]
    t70_celadonite[i] = bit[0,-3,24]

    


#############################
# SIMULATIONS WITH MINERALS #
# MIXING RATIO 60/40        #
#############################

# load temps
t02 = np.loadtxt('r60t02.txt')
t04 = np.loadtxt('r60t04.txt')
t06 = np.loadtxt('r60t06.txt')
t08 = np.loadtxt('r60t08.txt')
t10 = np.loadtxt('r60t10.txt')
t12 = np.loadtxt('r60t12.txt')
t14 = np.loadtxt('r60t14.txt')
t16 = np.loadtxt('r60t16.txt')
t18 = np.loadtxt('r60t18.txt')
t20 = np.loadtxt('r60t20.txt')

t = [[t02[:,:]], [t04[:,:]], [t06[:,:]], [t08[:,:]], [t10[:,:]],
     [t12[:,:]], [t14[:,:]], [t16[:,:]], [t18[:,:]], [t20[:,:]]]

# grab for each experiment
t60_alk = np.zeros((10))
t60_glass = np.zeros((10))
t60_kaolinite = np.zeros((10))
t60_stilbite = np.zeros((10))
t60_saponite = np.zeros((10))
t60_albite = np.zeros((10))
t60_celadonite = np.zeros((10))
for i in range(len(temps)):
    bit = np.asarray(t[i])
    t60_alk[i] = bit[0,-2,3] - .002*.4
    t60_glass[i] = bit[0,-2,57]
    t60_kaolinite[i] = bit[0,-3,18]
    t60_stilbite[i] = bit[0,-3,14]
    t60_saponite[i] = bit[0,-3,22]
    t60_albite[i] = bit[0,-3,20]
    t60_celadonite[i] = bit[0,-3,24]






fig=plt.figure()

plt.rc('xtick', labelsize=8) 
plt.rc('ytick', labelsize=8)



#########################
# PLOT ALK FOR MINERALS #
#########################

ax = plt.subplot(2,2,1)


p = plt.plot(temps,t90_alk,'k^-',linewidth=2,label='mixing 90/10')
p = plt.plot(temps,t80_alk,'r^-',linewidth=2,label='mixing 80/20')
p = plt.plot(temps,t70_alk,'y^-',linewidth=2,label='mixing 70/30')
p = plt.plot(temps,t60_alk,'c^-',linewidth=2,label='mixing 60/40')

plt.title('ALKALINITY',fontsize=8)
plt.ylabel('ALK TO OCEAN [mol kgw$^{-1}$ yr$^{-1}$ per 1000 mol basalt]',
           fontsize=6)
plt.xlabel('T [$^{\circ}$C]',fontsize=10)

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.legend(handles, labels,loc='best',prop={'size':6}, ncol=1)



###########################
# PLOT GLASS FOR MINERALS #
###########################

ax = plt.subplot(2,2,2)

p = plt.plot(temps,t90_glass,'k^-',linewidth=2,label='mixing 90/10')
#p = plt.plot(temps,t60_glass,'r^-',linewidth=2,label='mixing 80/20')


plt.title('REMAINING BASALT',fontsize=8)
plt.ylabel('ALKALINITY',fontsize=6)
plt.xlabel('T [$^{\circ}$C]',fontsize=6)

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.legend(handles, labels,loc='best',prop={'size':6}, ncol=1)




########################
# PLOT ALK NO MINERALS #
########################

ax = plt.subplot(2,2,4)


plt.title('',fontsize=10)
plt.ylabel('',fontsize=6)
plt.xlabel('T [$^{\circ}$C]',fontsize=10)

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.legend(handles, labels,loc='best',prop={'size':6}, ncol=1)



##############
# PLOT CLAYS #
##############

ax = plt.subplot(2,2,3)

p = plt.plot(temps,t90_celadonite,'k-',linewidth=1,label='mixing 90/10')
p = plt.plot(temps,t90_stilbite,'k--',linewidth=1)
p = plt.plot(temps,t90_kaolinite,'k:',linewidth=2)
p = plt.plot(temps,t90_albite,'k-.',linewidth=1)
p = plt.plot(temps,t90_saponite,'k-',linewidth=2)


##p = plt.plot(temps,t60_celadonite,'c-',linewidth=1,label='mixing 60/40')
##p = plt.plot(temps,t60_stilbite,'c--',linewidth=1)
##p = plt.plot(temps,t60_kaolinite,'c:',linewidth=2)
##p = plt.plot(temps,t60_albite,'c-.',linewidth=1)
##p = plt.plot(temps,t60_saponite,'c-',linewidth=2)


#p = plt.plot(temps,t80_celadonite,'r-',linewidth=3,label='mixing 80/20')
#p = plt.plot(temps,t80_stilbite,'r--',linewidth=3)


plt.ylabel('[mol / 2kyr]',fontsize=6)
plt.xlabel('T [$^{\circ}$C]',fontsize=10)

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.legend(handles, labels,loc='best',prop={'size':6}, ncol=1)



###########################
# PLOT 2,2 #
###########################

ax = plt.subplot(2,2,4)



plt.title('',fontsize=8)
plt.ylabel('',fontsize=6)
plt.xlabel('',fontsize=6)

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])
plt.legend(handles, labels,loc='best',prop={'size':8}, ncol=2)


plt.subplots_adjust(hspace=.25, wspace=.25)
plt.savefig('flushExp.png')
